import io
import os
import re

import numpy
import numpy as np
import pandas
import pandas as pd
from Bio import AlignIO, SeqIO

import settingmain


class Sequence:
    def __init__(self, seq, metadata=None):
        self.seq = seq.upper()
        if metadata:
            self.metadata = metadata
        else:
            self.metadata = dict()

    def find_with_regex(self, motif, ignore=None):
        pattern = re.compile(motif)
        new_str = ""
        if ignore is not None:
            for i in range(len(ignore)):
                if not ignore[i]:
                    new_str += self.seq[i]

        for i in pattern.finditer(new_str):
            for m in range(1, len(i.groups()) + 1):
                yield slice(i.start(m), i.end(m))

    def gaps(self):
        s = np.arange(len(self.seq))
        for i in range(len(s)):
            if self.seq[i] == '-':
                s[i] = True
            else:
                s[i] = False
        return s

    def count(self, char, start, end):
        return self.seq.count(char, start, end)

    def __repr__(self):
        return self.seq

    def __str__(self):
        return self.seq


class ProteinCollection:
    regex_motif = '(?=(N[^PX][ST]))'

    def __init__(self, alignment_fn=None, multifasta_fn=None, keeping=None):
        """

        :type multifasta_fn: str
        Multifasta filename
        :type alignment_fn: str
        Alignment filename
        """
        self.alignment_fn = alignment_fn
        self.multifasta_fn = multifasta_fn
        self.keeping = keeping
        self.motifs = []
        self.seqs = {}
        self.parse_motifs()

    def parse_motifs(self):
        if self.alignment_fn:

            alignment = dict()

            if type(self.alignment_fn) == str:
                align_gen = AlignIO.read(self.alignment_fn, 'clustal')
            else:
                align_gen = AlignIO.read(io.BytesIO(self.alignment_fn), 'clustal')
            for rec in align_gen:
                p = Sequence(str(rec.seq).upper(), metadata={'id': rec.id})
                p.metadata['motifs'] = [match for match in p.find_with_regex(self.regex_motif, ignore=p.gaps())]
                p.metadata['motif_fin'] = {}
                for i in p.metadata['motifs']:
                    gap_count = p.count('-', start=0, end=i.start)
                    p.metadata['motif_fin'][i.start - gap_count] = gap_count
                alignment[rec.id] = p
            self.motifs = alignment

            self.seqs = {rec.id: Sequence(str(rec.seq).upper(), metadata={'id': rec.id}) for rec in
                         SeqIO.parse(self.multifasta_fn, 'fasta')}
            print(self.seqs)

        else:
            seq = dict()
            for rec in SeqIO.parse(self.multifasta_fn, 'fasta'):
                #print(rec.id)
                if self.keeping:
                    if rec.id in self.keeping:
                        pass
                    else:
                        continue
                p = Sequence(str(rec.seq).upper(), metadata={'id': rec.id})

                p.metadata['motifs'] = [match for match in p.find_with_regex(self.regex_motif, ignore=p.gaps())]
                p.metadata['motif_fin'] = {}
                for i in p.metadata['motifs']:
                    gap_count = p.count('-', 0, i.start)
                    p.metadata['motif_fin'][i.start - gap_count] = gap_count
                seq[rec.id] = p
                self.seqs[rec.id] = p
            self.motifs = seq

    def run(self):
        pass


class Analyses:
    _mods = ['HexNAc', 'HexNAc(1)dHex(1)', 'Oxidation', 'Deamidated', 'N-Term(Gln->pyro-Glu)']

    def __init__(self, excel_fn, output_fn, protein_n, protein_collection, pattern, zip_handler=None, job_id=None):
        """

        :type protein_collection: ProteinCollection
        :type output_fn: str
        :type excel_fn: str

        """

        self.output_fn = output_fn
        self.protein_n = protein_n
        self.protein = protein_collection.seqs[protein_n]
        self.alignment = protein_collection.motifs[protein_n]
        # print(self.alignment.metadata)
        self.regex_motif = pattern
        self.zip_handler = zip_handler
        self.job_id = job_id
        if self.job_id is not None:
            self.output_folder = os.path.join(settingmain.APP_TEMP, job_id)
        else:
            self.output_folder = settingmain.APP_TEMP
        # print(self.df)
        if type(excel_fn) == str:
            self.df = pd.read_excel(excel_fn)
        else:
            self.df = excel_fn
        self.origin_df_columns_length = len(self.df.columns)

    def process(self):

        columns = ('aa number', 'Expanded aa', 'Expanded motif', 'Number of N-linked sequons', 'Number of Asn residues',
                   'Total number of H modification', 'Total number of Asn modification',
                   'Total unmodified Asn Residues')
        arr = []
        ind = ()
        row = 0
        # Select appropriate columns for loop of the protein of interset

        for _, r in self.df.iterrows():
            seq = r['Annotated Sequence']
            mod = r['Modifications']
            # print(i3)

            p = Sequence(seq.upper())
            ind += (row,)
            # print(p)
            # Finding N-glycosylation motifs using regex in the Annotated Sequence
            p.metadata['motifs'] = [match for match in p.find_with_regex(self.regex_motif, ignore=p.gaps())]
            # Finding the position of the fragment on the original protein sequence
            origin = self.protein.seq.index(p.seq)
            # Mapping the positions of the position of each motifs identified to their positions on the original
            # protein sequence
            p.metadata['origin_map'] = [i.start + origin for i in p.metadata['motifs']]
            p.metadata['origin'] = str(origin + 1) + '-' + str(origin + len(p.seq))
            # Expand the end of the fragments through mapping to the original protein sequence
            p.metadata['expanded_window'] = Sequence(self.protein.seq[origin:origin + len(p.seq) + 2])
            # Identify motifs within the expanded fragments
            p.metadata['expanded_window_motifs'] = [match for match in
                                                    p.metadata['expanded_window'].find_with_regex(self.regex_motif,
                                                                                                  ignore=p.metadata[
                                                                                                      'expanded_window']
                                                                                                  .gaps())]
            p.metadata['origin_map'] = [i.start + origin for i in p.metadata['expanded_window_motifs']]
            p.metadata['origin'] = str(origin + 1) + '-' + str(origin + len(p.seq))
            # print(result_dict)
            # print(self.alignment.metadata['motif_fin'])
            # print(p.metadata['origin_map'])

            result_dict = {
                'aa number': p.metadata['origin'],
                'Number of N-linked sequons': len(p.metadata['expanded_window_motifs']),
                'Expanded aa': str(p.metadata['expanded_window'].seq[-2:])
            }
            if len(p.metadata['expanded_window_motifs']) > len(p.metadata['motifs']):
                result_dict.update({'Expanded motif': str(
                    p.metadata['expanded_window'].seq[p.metadata['expanded_window_motifs'][-1]])})
            result_dict.update({
                'Number of Asn residues': p.count('N', 0, len(p.seq)),
            })

            s = []
            count_HexNAc = 0
            count_deamin = 0
            asn_dict = dict()
            if not pd.isnull(mod):
                s = mod.split(';')

                for mod in s:
                    rpattern = re.compile('(N(\w+))')
                    npos = rpattern.search(mod)
                    if npos:
                        # if 'HexNAc' in mod:
                        if 'Hex' in mod:
                            count_HexNAc += 1
                            asn_dict[str(npos.group(2))] = 'H'
                        elif 'Deamidated' in mod:
                            count_deamin += 1
                            asn_dict[str(npos.group(2))] = 'D'

            result_dict.update({
                'Total number of H modification': count_HexNAc,
                'Total number of Asn modification': count_HexNAc + count_deamin,
                'Total unmodified Asn Residues': p.count('N', 0, len(p.seq)) - (count_HexNAc + count_deamin)
            })

            # Matching of sequon characteristics
            n_sequon = 1
            extra = 1
            for n in p.metadata['origin_map']:
                if n in self.alignment.metadata['motif_fin']:
                    if '{0}. -linked alignment'.format(str(n_sequon)) not in columns:
                        columns += ('{0}. -linked alignment'.format(str(n_sequon)),)
                    if '{0}. -linked protein'.format(str(n_sequon)) not in columns:
                        columns += ('{0}. -linked protein'.format(str(n_sequon)),)
                    if '{0}. Match'.format(str(n_sequon)) not in columns:
                        columns += ('{0}. Match'.format(str(n_sequon)),)
                    result_dict.update({
                        '{0}. -linked alignment'.format(str(n_sequon)): 'N{0}'.format(
                            str(n + self.alignment.metadata['motif_fin'][n] + 1)),
                        '{0}. -linked protein'.format(str(n_sequon)): 'N{0}'.format(
                            str(n + 1)),
                        # '{0}. Match'.format(str(n_sequon)):
                    })
                    if result_dict['Number of N-linked sequons'] == count_HexNAc:
                        result_dict.update({
                            '{0}. Match'.format(str(n_sequon)): 'H'
                        })
                    else:

                        if str(n - origin + 1) in asn_dict:
                            if asn_dict[str(n - origin + 1)] == 'D':
                                if not count_HexNAc:
                                    result_dict.update({
                                        '{0}. Match'.format(str(n_sequon)): 'D'
                                    })
                                else:
                                    result_dict.update({
                                        '{0}. Match'.format(str(n_sequon)): 'D/H'
                                    })
                                    if result_dict['Total unmodified Asn Residues']:
                                        result_dict['{0}. Match'.format(str(n_sequon))] += '/U'
                            else:
                                if not count_deamin:
                                    result_dict.update({
                                        '{0}. Match'.format(str(n_sequon)): 'H'
                                    })
                                else:
                                    result_dict.update({
                                        '{0}. Match'.format(str(n_sequon)): 'D/H'
                                    })
                                if result_dict['Total unmodified Asn Residues']:
                                    result_dict['{0}. Match'.format(str(n_sequon))] += '/U'
                        else:
                            if count_deamin and count_HexNAc:
                                result_dict.update({
                                    '{0}. Match'.format(str(n_sequon)): 'D/H'
                                })
                                if result_dict['Total unmodified Asn Residues']:
                                    result_dict['{0}. Match'.format(str(n_sequon))] += '/U'
                            elif count_HexNAc:
                                result_dict.update({
                                    '{0}. Match'.format(str(n_sequon)): 'H'
                                })
                                if result_dict['Total unmodified Asn Residues'] and (
                                        len(p.metadata['origin_map']) > result_dict['Total number of H modification']):
                                    result_dict['{0}. Match'.format(str(n_sequon))] += '/U'
                            else:
                                result_dict.update({
                                    '{0}. Match'.format(str(n_sequon)): 'U'
                                })

                n_sequon += 1
                # print(result_dict)
            arr.append(result_dict)
            # print(result_dict)
            row += 1

        extend_df = pd.DataFrame(arr, columns=columns, index=self.df.index)
        self.processed_df = pd.concat([self.df, extend_df], axis=1)

    def compile(self, sites_number: int, separate_h: bool = False, glycans=None, aggregate=None) -> object:
        temp = self.processed_df.copy()
        length = len(temp.columns)
        temp['Annotated Sequence'] = [x.upper() for x in temp['Annotated Sequence']]
        aggregate_label = ""
        if aggregate:
            if len(aggregate) > 0:
                aggregate_label = "-".join(aggregate)
                for r in range(self.origin_df_columns_length + 8, length, 3):
                    temp[temp.columns[r + 2]].replace(aggregate, aggregate_label, inplace=True)
                    # print(temp[temp.columns[r+2]].replace(aggregate, aggregate_label))
        else:
            aggregate = []
        base = ['H', 'D', 'U']
        base_columns = ["Annotated Sequence", "Charge", "Modifications"]
        additional_columns = []
        if "XCorr" in temp.columns:
            additional_columns.append("Area")
        if separate_h:
            additional_columns.append("Glycan composition")
        gr = [i for i in base_columns + additional_columns if i in temp.columns]
        position = dict()
        unique = []
        glycan_comp = []

        temp["Modifications"] = temp["Modifications"].fillna("None")

        if separate_h:
            if "Glycan composition" in gr:
                temp["Glycan composition"] = temp["Glycan composition"].fillna("None")
                glycan_comp = temp["Glycan composition"].unique()
                if len(glycan_comp) > 0:
                    if not glycans:
                        glycans = base
                    if len(glycans) == 0:
                        glycan_comp = np.append(glycan_comp,
                                                [i for i in ['D', 'U'] if i not in aggregate])

                    else:
                        glycan_comp = np.append(glycan_comp,
                                                [i for i in base if i in glycans if i != 'H' if i not in aggregate])

                else:
                    glycan_comp = [i for i in base if i not in aggregate]
        else:
            glycan_comp = [i for i in base if i not in aggregate]

        if aggregate:
            glycan_comp = np.append(glycan_comp, [aggregate_label])

        for i in temp.groupby(gr):
            if 'XCorr' in temp.columns:
                unique_row = i[1].loc[i[1]['XCorr'].idxmax()]
            else:
                unique_row = i[1].loc[i[1]['Area'].idxmax()]
            if (unique_row['Number of N-linked sequons']) <= sites_number:
                for r in range(self.origin_df_columns_length + 8, length, 3):
                    if pd.notnull(unique_row[r]):
                        if unique_row[r + 2] in glycans:
                            if unique_row[r] not in position:
                                position[unique_row[r]] = dict()
                            if unique_row[r + 1] not in position[unique_row[r]]:
                                position[unique_row[r]][unique_row[r + 1]] = dict(total=0)

                            if not separate_h:

                                if unique_row[r + 2] not in position[unique_row[r]][unique_row[r + 1]]:
                                    position[unique_row[r]][unique_row[r + 1]][unique_row[r + 2]] = 0
                                position[unique_row[r]][unique_row[r + 1]][unique_row[r + 2]] += unique_row['Area']
                                position[unique_row[r]][unique_row[r + 1]]['total'] += unique_row['Area']
                            else:
                                if "Glycan composition" in gr and separate_h:
                                    if 'H' in unique_row[r + 2]:
                                        if unique_row["Glycan composition"] not in position[unique_row[r]][
                                            unique_row[r + 1]]:
                                            position[unique_row[r]][unique_row[r + 1]][unique_row["Glycan composition"]] = 0
                                        position[unique_row[r]][unique_row[r + 1]][unique_row["Glycan composition"]] += \
                                            unique_row['Area']
                                        position[unique_row[r]][unique_row[r + 1]]['total'] += unique_row['Area']
                                    else:
                                        if unique_row[r + 2] not in position[unique_row[r]][unique_row[r + 1]]:
                                            position[unique_row[r]][unique_row[r + 1]][unique_row[r + 2]] = 0
                                        position[unique_row[r]][unique_row[r + 1]][unique_row[r + 2]] += unique_row['Area']
                                        position[unique_row[r]][unique_row[r + 1]]['total'] += unique_row['Area']
                                else:
                                    if unique_row[r + 2] not in position[unique_row[r]][unique_row[r + 1]]:
                                        position[unique_row[r]][unique_row[r + 1]][unique_row[r + 2]] = 0
                                    position[unique_row[r]][unique_row[r + 1]][unique_row[r + 2]] += unique_row['Area']
                                    position[unique_row[r]][unique_row[r + 1]]['total'] += unique_row['Area']
                unique.append(unique_row)
        unique_df = pd.DataFrame(unique)
        self.save(os.path.join(self.output_folder, 'unique' + self.output_fn), df=unique_df)
        if self.zip_handler is not None:
            self.zip_handler.write(os.path.join(self.output_folder, 'unique' + self.output_fn),
                                   'unique' + self.output_fn)
        columns = ['position'] + [i for i in glycan_comp if pd.notnull(i) and i != "None"]
        summary = {i: [] for i in glycan_comp if pd.notnull(i) and i != "None"}
        # if len(glycan_comp) > 0:
        #
        #     summary = {i: [] for i in glycan_comp if pd.notnull(i) and i != "None"}
        # else:
        #
        #     summary = {'H': [], 'D': [], 'U': []}
        ind = []
        for i in position:
            for i2 in position[i]:
                ind.append(i + "(" + i2 + ")")
                for s in summary:
                    if s in position[i][i2]:
                        summary[s].append(position[i][i2][s] / position[i][i2]['total'])
                    else:
                        summary[s].append(np.nan)
        summary['position'] = ind
        print(position)
        summary = pd.DataFrame(summary, columns=columns)
        # self.save(os.path.join(settingmain.APP_TEMP, 'summary'+self.output_fn), df=summary)
        return summary

    def save(self, filename='output.xlsx', df=None, index=False):
        if df is None:
            df = self.df
        writer = pd.ExcelWriter(filename)
        df.to_excel(writer, 'Sheet1', index=index)
        writer.save()


def analyze(b, cat, dfs, job_id, p, wdf, zf):
    header = pd.MultiIndex.from_product([cat, b['reps']], names=["mod", "rep"])
    result = dict()
    for _, r in wdf.iterrows():
        d = Analyses(dfs[r['filename']], r['filename'], r['protein'], p, '(?=(N[^PX][ST]))', zip_handler=zf,
                     job_id=job_id)
        d.process()
        summary = d.compile(b['maxSites'])
        if r['protein'] + r['condition'] not in result:
            result[r['protein'] + r['condition']] = dict()

        if r['replicate'] not in result[r['protein'] + r['condition']]:
            result[r['protein'] + r['condition']][r['replicate']] = summary
        if 'pos' not in result[r['protein'] + r['condition']]:
            result[r['protein'] + r['condition']]['pos'] = summary['position'].values
        else:
            result[r['protein'] + r['condition']]['pos'] = np.concatenate(
                (result[r['protein'] + r['condition']]['pos'], summary['position'].values))
    return header, result


def analyze2(df_map, work_df, p, query, job_id, zf):
    summaries = []
    for _, r in work_df.iterrows():
        d = Analyses(df_map[r['filename']][r['protein']], r['filename'], r['protein'], p, '(?=(N[^PX][ST]))', zip_handler=zf,
                     job_id=job_id)
        d.process()
        summary = d.compile(query['maxSites'], query['separate_h'], glycans=query['glycans'],
                            aggregate=query['aggregation'])
        m_summary = summary.melt(['position'], summary.columns[1:])
        m_summary["reps"] = pd.Series([r["replicate"]] * len(m_summary.index), index=m_summary.index)
        m_summary["condition"] = pd.Series([r["condition"] + "_" + r["protein"]] * len(m_summary.index),
                                           index=m_summary.index)
        summaries.append(m_summary)
    results = pd.concat(summaries, ignore_index=True)
    results['glycan_composition'] = results['variable']
    results['experiments'] = results['value']

    return results
