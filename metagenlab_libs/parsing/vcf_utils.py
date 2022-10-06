#!/usr/bin/env python
import pandas
import vcf

class VCF():

    def __init__(self, vcf_file, info_prefix='', sample_id=None):

        self.vcf_reader =  vcf.Reader(filename=vcf_file)
        self.calls = [i for i in self.vcf_reader]
        if len(self.calls) > 0:
            self.info_df = self.parse_info(info_prefix, sample_id)
        else:
            self.info_df = None 
        self.info_id2description = {y.id: y.desc for x,y in self.vcf_reader.infos.items()}


    def parse_info(self, prefix='', sample_id=None):

        header = ["var"] + [f'{prefix}{i}' for i in self.calls[0].INFO]

        rows = []    
        for call in self.calls:
            var_name = f'{call.REF}{call.POS}{call.ALT[0]}'
            info_data = [call.INFO[i][0] if isinstance(call.INFO[i], list) else call.INFO[i] for i in call.INFO] 
            rows.append([var_name] + info_data)

        df = pandas.DataFrame(rows, columns=header)

        if sample_id:
            df["sample_id"] = [sample_id]*len(df)
            return df.set_index(["sample_id","var"])
        else:
            return df.set_index(["var"])

