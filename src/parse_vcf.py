
import pandas as pd
import numpy as np
import vcf

import warnings
import os
import sys

from .utility import (infile_handler,
                     outfile_handler,
                     is_gzip)



#################################
# Class defintion : Dataset
#################################
class Dataset(object):
    def __init__(self, sample_ids):
        self.sample_ids = list(sample_ids)
        self.data = {sample_id : Sample(sample_id) for sample_id in sample_ids}
        self.dbsnp_vcf = None


    def add_vcf(self, sample_id, vcfname, vcfpath):
        self.data[sample_id].add_vcf(vcfname, vcfpath)
        return

    def get_pass_variants(self):
        for sample_id in self.sample_ids:
            self.data[sample_id].get_pass_variants()
        return

    def get_consensus_union(self):
        for sample_id in self.sample_ids:
            pass
        pass


    def __str__(self):
        string = "Sample IDs : %s " % ', '.join(self.sample_ids)
        return string

    def __setitem__(self, sample_id, data):
        self.data[sample_id] = data
        return

    def __getitem__(self, sample_id):
        return self.data[sample_id]


#################################
# Class defintion : Sample
#################################
class Sample(object):
    def __init__(self, sample_id, metadata=None):
        self.sample_id = sample_id
        self.metadata = metadata
        self.vcfs = {}
        self.pass_vcfs = {}

    def process_sample(self):
        self.get_pass_variants()
        self.compute_consensus_variants()
        self.get_sanger_VAF()
        return

    def add_vcf(self, vcfname, vcfpath):
        self.vcfs[vcfname] = VCF(vcfpath, pandas_engine='c')
        return self.vcfs[vcfname]

    def get_pass_variants(self):
        filter_has_pass = lambda df : df.loc[map(lambda f: "PASS" in f,
                                                 df["FILTER"])]
        for vcfname in self.vcfs:
            self.pass_vcfs[vcfname] = filter_has_pass(self.vcfs[vcfname].df)
        return

    def get_vcf(self, vcfname):
        return self.vcfs[vcfname]

    def get_sample_type(self):
        return self.metadata["Sample Type"]

    def compute_consensus_variants(self):
        pass_vcfnames = sorted(self.pass_vcfs.keys())
        fst_name = pass_vcfnames[0]
        df_result = self.pass_vcfs[fst_name].loc[:, ["ALT"]]
        df_result.columns = ["ALT_%s" % fst_name]
        for name in pass_vcfnames[1:]:
            df2 = self.pass_vcfs[name].loc[:, ["ALT"]]
            df2.columns = ["ALT_%s" % name]
            df_result = pd.merge(df_result, df2,
                                 how="inner",
                                 on=["CHROM", "POS", "REF"])
        self.df_consensus = df_result
        return df_result

    def get_sanger_VAF(self, df_sanger=None, sanger_samplename="TUMOUR"):
        if not df_sanger:
            df_sanger = self.pass_vcfs["SANGER"]
        df_result = pd.DataFrame(index=df_sanger.index,
                                 columns=["DP", "VAF", "FVAF", "RVAF"])
        for index in df_sanger.index:
            row = df_sanger.loc[index, :]
            dp = VCF.get_row_info(row)["DP"]
            alt = row["ALT"]
            fmt = VCF.get_row_sample_genotype(row, sanger_samplename)
            fad = fmt['F' + alt + 'Z']
            rad = fmt['R' + alt + 'Z']
            ad = fad + rad
            df_result.loc[index, ["DP", "VAF", "FVAF", "RVAF"]] = [dp, ad / dp, fad / dp, rad / dp]
        self.df_sanger_VAF = df_result
        return df_result

    def __str__(self):
        string = "Sample ID : %s \n" % self.sample_id
        string += "VCFs     : %s \n" % str(self.vcfs.keys())
        string += "Metadata : \n%s \n" % str(self.metadata)
        return string

    def dereference_original_vcfs(self):
        self.vcfs = {}
        return

    def dereference_pass_vcfs(self):
        self.pass_vcfs = {}
        return


    def consensus_nondbSNP(self, dbsnp_snv_variant_dict):
        vcf_names = self.vcfs.keys()
        def in_dbSNP(row):
            alts = np.array(row[["ALT_%s" % name for name in vcf_names]])
            alt = alts[0]
            if np.all(alts == alt):
                var = row.name[0] + '_' + str(row.name[1]) + '_' + row.name[2] + '>' + alt
                if var in dbsnp_snv_variant_dict:
                    return True
            else:
                warnings.warn("WARNING: disagreeing consensus calls: %s\n%s"
                              %(self.sample_id, str(row)))
            return False
        df = self.df_consensus.copy(deep=True)
        df["in_dbSNP"] = df.apply(in_dbSNP, axis=1)
        self.df_consensus_nondbSNP_VAF = (
            pd.merge((df.groupby("in_dbSNP")
                        .get_group(False)
                        .drop("in_dbSNP", axis=1)),
                     self.df_sanger_VAF,
                     how="inner",
                     on=["CHROM", "POS", "REF"]))
        return








#################################
# Class defintion : VCF
#################################
class VCF(object):
    """
    This looks a lot like reimplementing PyVCF.Reader
    just with Pandas.DataFrame object to hold the main body
    of the vcf file.
    """
    def __init__(self, vcf_filepath,
                 pandas_engine="c"):
        self.vcf_filepath = vcf_filepath
        self.pandas_engine = pandas_engine
        # parse the body of vcf into a Pandas.DataFrame object
        (self.df,
         self.header,
         self.SAMPLES) = VCF.parse_vcf(vcf_filepath,
         pandas_engine=pandas_engine)



    @staticmethod
    def parse_vcf(vcf_filepath, pandas_engine="python"):
        with infile_handler(vcf_filepath) as vcf_file:
            header = []
            for line in vcf_file:
                line = line.rstrip()
                if line.startswith("##"):
                # if it is a header line
                    header.append(line)
                elif line.startswith("#CHROM"):
                # if it is the column names line
                    columns = line[1:].split('\t')
                    # immediately break the loop
                    # so that the filehandle starts with the main body
                    # of the vcf file
                    break
            # which we then read into a Pandas.DataFrame object
            df = pd.read_csv(vcf_file, sep='\t',
                             header=None,
                             engine=pandas_engine,
                             dtype={0 : 'category',
                                    1 : np.uint32})
        # Give correct header
        df.columns = columns
        # slice out samples
        samples = columns[9:]
        # index dataframe using CHROM and POS (i.e. genomic loci)
        df.set_index(["CHROM", "POS", "REF"], inplace=True)
        # sort dataframe using genomic loci.
        df.sort_index(inplace=True)
        return df, header, samples

    def get_format_tags_intersection(self):
        common = set(self.df["FORMAT"][0].split(':'))
        for tags in self.df["FORMAT"]:
            current = set(tags.split(':'))
            common = common.intersection(current)
        return common

    def get_info_tags_intersection(self):
        get_tags = lambda info : set(field.split('=', maxsplit=1)[0]
                                     for field in info.split(';'))
        common = set(get_tags(self.df["INFO"][0]))
        for info in self.df["INFO"]:
            current = get_tags(info)
            common = common.intersection(current)
        return common

    def melt_genotype_column(self, samples=None):
        if samples is None:
            samples = self.SAMPLES
        elif type(samples) == int:
            samples = [self.SAMPLES[samples]]

        common_format_tags = self.get_format_tags_intersection()
        for sample in samples:
            for tag in common_format_tags:
                data = []
                col_name = sample + '_' + "FORMAT" + '_' + tag
                self.df[col_name] = np.zeros(self.df.shape[0])
                for i, row in self.df.iterrows():
                    row_sample_genotype = VCF.get_row_sample_genotype(row,
                                                                      sample)
                    data.append(row_sample_genotype[tag])
                self.df[col_name] = data
        return self.df

    def melt_info_column(self):
        common_info_tags = self.get_info_tags_intersection()
        for tag in common_info_tags:
            data = []
            for i, row in self.df.iterrows():
                data.append(VCF.get_row_info(row)[tag])
            self.df["INFO" + '_' + tag] = data
        return self.df

    @staticmethod
    def get_all_allelic_change(df):
        data = []
        for i, row in df.iterrows():
            for alt in row["ALT"].split(','):
                change = "{ref}>{alt}".format(ref=row["REF"], alt=alt)
                record = [row["CHROM"], row["POS"], change]
                data.append(record)
        return pd.DataFrame(data, columns=["CHROM", "POS", "CHANGE"])

    ########
    # VCF dataframe utilities
    ########
    @staticmethod
    def get_row_sample_genotype(row, sample_name):
        tags = row["FORMAT"].split(':')
        tag_values = row[sample_name].split(':')
        result = dict(zip(tags, tag_values))
        for key, val in result.items():
            try:
                result[key] = np.float(val)
            except ValueError:
                pass
        return result

    @staticmethod
    def get_row_info(row):
        result = {}
        for field in row["INFO"].split(';'):
            if '=' in field:
                key, val = field.split('=')
                result[key] = val
            else:
                result[field] = True
        for key, val in result.items():
            if type(val) is bool:
                continue
            try:
                result[key] = np.float(val)
            except ValueError:
                pass
        return result




############# Legacy.....

class Sanger_VCF(VCF):
    def __init__(self, vcf_filepath, pandas_engine="c"):
        # Parent class initiation
        VCF.__init__(self, vcf_filepath, pandas_engine)

        self.vcf_reader = vcf.Reader(filename=vcf_filepath,
                                     compressed=is_gzip(vcf_filepath))
        if store_record:
            self.records = [rec for rec in self.vcf_reader]
        else:
            self.records = None

        self.FORMAT = self.vcf_reader.formats
        self.INFO = self.vcf_reader.infos
        self.SAMPLES = self.vcf_reader.samples








