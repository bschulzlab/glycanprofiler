import unittest
import pandas as pd
import os
from common import ProteinCollection, Analyses

p = ProteinCollection(multifasta_fn=r"C:\Users\localadmin\PycharmProjects\glycanProfiler\temp\Schulz_RP_WT_FLAG_HuNCCAM_20_718_8His_HuPST.FASTA")
filenames = ["ruby-1.xlsx", "ruby-2.xlsx", "ruby-3.xlsx"]
tempdf = {}
for i in filenames:
    df = pd.read_excel(os.path.join("temp", i))
    df = df[(df['Master Protein Accessions'].notnull()) &
                                df['Master Protein Accessions'].str.contains("RP2908_WT_FLAG_HuNCAM_20-718_8His") & (
                                        ~df['Area'].isnull()) & (df['Search Engine Rank'] == 1)]
    tempdf[i] = df


class TestAnalysesMethods(unittest.TestCase):

    def test_process(self):
        for df_name in tempdf:
            d = Analyses(tempdf[df_name], df_name, "RP2908_WT_FLAG_HuNCAM_20-718_8His", p, '(?=(N[^PX][ST]))',
                         job_id="1")
            d.process()
            self.assertGreater(len(d.processed_df.columns), len(tempdf[df_name].columns))

    def test_compile(self):
        for df_name in tempdf:
            d = Analyses(tempdf[df_name], df_name, "RP2908_WT_FLAG_HuNCAM_20-718_8His", p, '(?=(N[^PX][ST]))',
                         job_id="1")
            d.process()
            summary = d.compile(1)
            self.assertFalse(summary.empty)

    def test_analyze(self):
        summaries = []
        for df_name in tempdf:
            d = Analyses(tempdf[df_name], df_name, "RP2908_WT_FLAG_HuNCAM_20-718_8His", p, '(?=(N[^PX][ST]))',
                         job_id="1")
            d.process()
            summary = d.compile(1, True)
            m_summary = summary.melt(['position'], summary.columns[1:])
            m_summary["reps"] = pd.Series([df_name] * len(m_summary.index), index=m_summary.index)
            m_summary["condition"] = pd.Series(["a"] * len(m_summary.index), index=m_summary.index)
            summaries.append(m_summary)
        results = pd.concat(summaries, ignore_index=True)
        results['glycan_composition'] = results['variable']
        results['experiments'] = results['value']
        print(results)
        for i, d in results.groupby(["condition"]):
            fin = pd.pivot_table(d,
                                 values=['experiments'],
                                 index=['position'],
                                 columns=['glycan_composition', 'reps'],
                                 dropna=False)
            writer = pd.ExcelWriter('ruby.xlsx')
            fin.to_excel(writer, 'Sheet1')
            writer.save()

    def test_analyze_aggregate(self):
        summaries = []
        for df_name in tempdf:
            d = Analyses(tempdf[df_name], df_name, "RP2908_WT_FLAG_HuNCAM_20-718_8His", p, '(?=(N[^PX][ST]))',
                         job_id="1")
            d.process()
            summary = d.compile(1, True, aggregate=['D', 'U'])
            m_summary = summary.melt(['position'], summary.columns[1:])
            m_summary["reps"] = pd.Series([df_name] * len(m_summary.index), index=m_summary.index)
            m_summary["condition"] = pd.Series(["a"] * len(m_summary.index), index=m_summary.index)
            summaries.append(m_summary)
        results = pd.concat(summaries, ignore_index=True)
        results['glycan_composition'] = results['variable']
        results['experiments'] = results['value']
        print(results)
        for i, d in results.groupby(["condition"]):
            fin = pd.pivot_table(d,
                                 values=['experiments'],
                                 index=['position'],
                                 columns=['glycan_composition', 'reps'],
                                 dropna=False)
            writer = pd.ExcelWriter('ruby.xlsx')
            fin.to_excel(writer, 'Sheet1')
            writer.save()

    def test_analyze_no_comp(self):
        summaries = []
        for df_name in tempdf:
            d = Analyses(tempdf[df_name], df_name, "RP2908_WT_FLAG_HuNCAM_20-718_8His", p, '(?=(N[^PX][ST]))',
                         job_id="1")
            d.process()
            summary = d.compile(1, False, aggregate=['D', 'U'])
            m_summary = summary.melt(['position'], summary.columns[1:])
            m_summary["reps"] = pd.Series([df_name] * len(m_summary.index), index=m_summary.index)
            m_summary["condition"] = pd.Series(["a"] * len(m_summary.index), index=m_summary.index)
            summaries.append(m_summary)
        results = pd.concat(summaries, ignore_index=True)
        results['glycan_composition'] = results['variable']
        results['experiments'] = results['value']
        print(results)
        for i, d in results.groupby(["condition"]):
            fin = pd.pivot_table(d,
                                 values=['experiments'],
                                 index=['position'],
                                 columns=['glycan_composition', 'reps'],
                                 dropna=False)
            writer = pd.ExcelWriter('ruby.xlsx')
            fin.to_excel(writer, 'Sheet1')
            writer.save()