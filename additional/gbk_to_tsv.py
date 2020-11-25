#!/usr/bin/env python3

from ParseGBK import ParseGBK

"""
DESCRIPTION: This script extracts genes, locus_tags, protein products, E.C. numbers, and translation sequences from an
             S. auerus strain used to test a new iron-containing drug called FS-1.
"""


def main():
    gbk = ParseGBK('control-strain_NC.gbf')
    gbk.to_tsv(['gene', 'locus_tag', 'product', 'EC_number', 'translation'], 'CDS', 'control-strain_NC.tsv')


if __name__ == '__main__':
    main()
