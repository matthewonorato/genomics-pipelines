#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd


class ParseGBK:
    """
    Parses a single .gbk file using Biopython's SeqIO.read function --> only for GBK file with a single record
    """

    def __init__(self, gbk_file):
        """
        Instantiates a ParseGBK class object
        Usage (with import ParseGBK): gbk = ParseGBK.ParseGBK(gbk_file)
        :param gbk_file: single-record GenBank file
        """
        self.gbk_file = gbk_file

    def __open_gbk(self, gbk_file):
        """
        Opens a single-record GBK file and returns SeqRecord object (for use in other class methods)
        :param gbk_file: same as gbk_file in __init__
        :return: SeqRecord object
        """
        with open(gbk_file, 'r') as handle:
            return SeqIO.read(handle, 'genbank')

    def __str__(self):
        """
        Displays basic info about genome, such as name, accession number, and strain
        :return: string of record
        """
        record = self.__open_gbk(self.gbk_file)
        return str(record)

    def __len__(self):
        """
        Displays entire length of genome
        :return: integer value
        """
        record = self.__open_gbk(self.gbk_file)
        return len(record)

    def sequence(self):
        """
        Displays DNA sequence as a gigantic string; possibly useful for strain-to-strain comparison of DNA differences
        :return: string of entire sequence
        """
        record = self.__open_gbk(self.gbk_file)
        return record.seq

    def features(self):
        """
        Displays all unique feature types found in GBK file
        :return: set of feature types
        """
        record = self.__open_gbk(self.gbk_file)
        return set([feature.type for feature in record.features])

    def get_feature(self, feat_type):
        """
        Extracts all feature values for a given feature type
        :param feat_type: feature type, such as 'CDS'
        :return: list of all feature values
                 (note: each element in the list is a SeqFeature object --> Ex: f[0].location.strand)
        """
        record = self.__open_gbk(self.gbk_file)
        return [feature for feature in record.features if feature.type == feat_type]

    def feature_locations(self, feat_type):
        """
        Extracts start & stop position and strand for all features of a given type
        :param feat_type: feature type, such as 'CDS'
        :return: list of maps; each map contains 'start', 'stop', and 'strand' keys for easy parsing
        """
        feature = self.get_feature(feat_type)
        return [{'start': f.location.start, 'stop': f.location.end, 'strand': f.strand} for f in feature]

    def qualifiers(self):
        """
        Displays all unique qualifiers found in GBK file
        :return: set of qualifiers (note: each feature contains it's own qualifiers)
        """
        record = self.__open_gbk(self.gbk_file)

        qualifiers = set()
        for feature in record.features:
            for qualifier in feature.qualifiers.keys():
                qualifiers.add(qualifier)
        return qualifiers

    def get_qualifier(self, qual_type):
        """
        Extracts all qualifier values for a given qualifier type
        :param qual_type: qualifier type, such as 'gene'
        :return: list of all qualifier values
        """
        record = self.__open_gbk(self.gbk_file)

        qual_vals = []
        for feature in record.features:
            for k, v in feature.qualifiers.items():
                if k == qual_type:
                    qual_vals += v
        return qual_vals

    def qualifier_of_feature(self, qual_type, feat_type):
        """
        Extracts a qualifier type from a feature (ex. genes in CDSs)
        :param qual_type: qualifier type, such as 'gene'
        :param feat_type: feature type, such as 'CDS'
        :return: list of qualifier type specific to a feature
        """
        self.__open_gbk(self.gbk_file)
        feature = self.get_feature(feat_type)

        q_of_f_vals = []
        for f in feature:
            if qual_type in f.qualifiers.keys():
                qual_vals = []
                for k, v in f.qualifiers.items():
                    if k == qual_type:
                        qual_vals += v
                if len(qual_vals) == 1:
                    qual_vals = qual_vals[0]
                q_of_f_vals.append(qual_vals)
            else:
                q_of_f_vals += ['-']

        return q_of_f_vals

    def to_df(self, quals, feat_type):
        """
        Creates dataframe of specified features from GBK file
        :param quals: list of qualifier types, such as ['gene', 'locus_tag', 'product', 'EC_number', 'translation']
        :param feat_type: feature type, such as 'CDS'
        :return: dataframe of extracted qualifiers
        """
        q_to_vals = {}

        for qual in quals:
            q_to_vals[qual] = self.qualifier_of_feature(qual, feat_type)
        df = pd.DataFrame(q_to_vals)
        print(df.describe())

        return df

    def to_tsv(self, quals, feat_type, tsv_path):
        """
        Creates tab-separated file (.tsv) of specified features from GBK file
        :param quals: list of qualifier types, such as ['gene', 'locus_tag', 'product', 'EC_number', 'translation']
        :param feat_type: feature type, such as 'CDS'
        :param tsv_path: path for written .tsv file
        :return: written .tsv file
        """
        df = self.to_df(quals, feat_type)
        df.to_csv(tsv_path, index=False, sep='\t', header=True, mode='w', encoding=None)
        return
