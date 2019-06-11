from tornado import web, gen, ioloop, escape
import pandas as pd
import io
import os
import numpy as np
import settingmain
import uuid
import zipfile

from common import ProteinCollection, analyze2


class BaseHandler(web.RequestHandler):
    def set_default_headers(self):
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept")
        self.set_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS, PUT')

    def options(self):
        self.set_status(204)
        self.finish()


class GlycanProfilerHandler(BaseHandler):
    def data_received(self, chunk):
        pass

    @gen.coroutine
    def post(self, *args, **kwargs):
        print(self.request)
        dfs = dict()
        work_array = []
        job_id = uuid.uuid4().hex
        serve_file = os.path.join(settingmain.APP_STATIC, job_id + '.zip')
        folder = os.path.join(settingmain.APP_TEMP, job_id)
        os.mkdir(folder)
        boil = escape.json_decode(self.request.files['boilerplate'][0].body)

        if 'alignment' in self.request.files:
            align = io.StringIO(self.request.files['alignment'][0].body.decode('UTF-8'))
        else:
            align = None
        mf = os.path.join(folder, 'multifasta.fasta')
        with open(mf, 'wb') as multifasta:
            multifasta.write(self.request.files['fasta'][0].body)
        # mf = io.StringIO(self.request.files['fasta'][0].body.decode('UTF-8'))
        p = ProteinCollection(alignment_fn=align, multifasta_fn=mf)
        for i in self.request.files['files']:
            if i.filename not in dfs:
                if i.filename.endswith('.xlsx'):
                    for b in boil:
                        if i.filename in b['repsMap']:
                            df = pd.read_excel(io.BytesIO(i.body))
                            print("Pre-processing: %s" % i.filename)
                            df = df[(df['Master Protein Accessions'].notnull()) & (
                                df['Master Protein Accessions'].str.contains(b['protein'])) & (
                                        ~df['Area'].isnull()) & (df['Search Engine Rank'] == 1) & (
                                            df['Area'] >= b['minimumArea'])]
                            if len(df.index) > 0:
                                dfs[i.filename] = df
                                work_array.append(
                                    [b['repsMap'][i.filename], b['condsMap'][i.filename], b['protein'], i.filename])
                                break
        work_df = pd.DataFrame(work_array, columns=['replicate', 'condition', 'protein', 'filename'])
        print(p)
        print(work_df)
        output = []
        for b in boil:
            wdf = work_df[work_df['protein'] == b['protein']]
            # cat = ["H", "D", "U"]
            # lrep = len(b['reps'])
            # lcat = len(cat)
            with zipfile.ZipFile(serve_file, 'w') as zf:
                results = analyze2(dfs, wdf, p, b, job_id, zf)
                for i, d in results.groupby(["condition"]):
                    fin = pd.pivot_table(d,
                                         values=['experiments'],
                                         index=['position'],
                                         columns=['glycan_composition', 'reps'],
                                         dropna=False)

                    writer = pd.ExcelWriter(os.path.join(folder, i + '.xlsx'))
                    fin.to_excel(writer, 'Sheet1')
                    writer.save()
                    zf.write(os.path.join(folder, i + '.xlsx'), i + '.xlsx')
                # for i in result:
                #     data = []
                #     result[i]['pos'] = sorted(np.unique(result[i]['pos']))
                #
                #     for i2 in result[i]['pos']:
                #         row = lrep * lcat * [np.nan]
                #         for lr in range(lrep):
                #             if b['reps'][lr] in result[i]:
                #                 print(b['reps'][lr], i2)
                #                 print(result[i][b['reps'][lr]]['position'])
                #                 print(result[i][b['reps'][lr]])
                #                 # if result[i][b['reps'][lr]]['position'].values != []:
                #                 if not result[i][b['reps'][lr]].empty:
                #                     new_sel = result[i][b['reps'][lr]][result[i][b['reps'][lr]]['position'] == i2]
                #                     for n in new_sel.iterrows():
                #                         for lc in range(lcat):
                #                             row[lr + lrep * lc] = n[1][cat[lc]]
                #         data.append(row)
                #
                #     fin_df = pd.DataFrame(data, columns=header, index=result[i]['pos'])
                #     writer = pd.ExcelWriter(os.path.join(folder, i + '.xlsx'))
                #     fin_df.to_excel(writer, 'Sheet1')
                #     writer.save()
                #     zf.write(os.path.join(folder, i + '.xlsx'), i + '.xlsx')

            output.append(dict(protein=[b['protein']], filename=job_id + '.zip'))
        print("Finished.")
        self.write(dict(data=output))


settings = {
    "autoreload": True,
    "debug": True,
}

if __name__ == "__main__":
    application = web.Application([
        (r"/api/gpp/upload/", GlycanProfilerHandler),
        (r"/static/(.*)", web.StaticFileHandler, {"path": "./static"}),
    ], **settings)
    application.listen(9001)
    ioloop.IOLoop.current().start()
