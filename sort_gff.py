#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "v1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_csv(file, sep='\t'):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def split_attr(attributes):

    r = collections.OrderedDict()
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if "=" not in content:
            print("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def get_id(attributes):

    gene_attr = split_attr(attributes)

    if 'ID' in gene_attr:
        gene_id = gene_attr['ID']
    elif 'Parent' in gene_attr:
        gene_id = gene_attr['Parent']
    else:
        raise Exception('Gene has no name.')

    return gene_id


def read_gff(file):

    gene_dict = {}
    mrna_dict = {}
    cds_dict = {}


    for line in read_csv(file):
        gene_id = get_id(line[-1])

        if line[2]=="exon":
            continue
        if line[2]=="mRNA":
            mrna_dict[gene_id] = line
        elif line[2]=='CDS':
            if gene_id not in cds_dict:
                cds_dict[gene_id] = []
            cds_dict[gene_id].append(line)
        else:
            gene_id = gene_id.strip('-gene')
            gene_dict[gene_id] = line

    return gene_dict, mrna_dict, cds_dict


def update_gene_name(gene_line, old_name, locustag, number):

    if locustag!='':
        new_name = '%s_%05d' % (locustag, number)
        gene_line[-1] = gene_line[-1].replace(old_name, new_name)

    return '\t'.join(gene_line)
    

def sort_gff(file, locustag):

    n = 0
    gene_dict, mrna_dict, cds_dict = read_gff(file)

    for line in sorted(gene_dict.items(), key = lambda x:(x[1][0], int(x[1][3]), int(x[1][4]))):
        n += 1
        print(update_gene_name(line[1], line[0], locustag, n))

        if line[0] in mrna_dict:
            print(update_gene_name(mrna_dict[line[0]], line[0], locustag, n))

        if line[0] not in  cds_dict:
            continue
        for cdsline in sorted(cds_dict[line[0]], key = lambda x: (x[0], int(x[3]), int(x[4]))):
            cds_line = update_gene_name(cdsline, line[0], locustag, n)
            cdsline[2] = "exon"
            print(update_gene_name(cdsline, line[0], locustag, n))
            print(cds_line)

    return 0


def add_help_args(parser):

    parser.add_argument('gff',
        help='Input genome annotation gff file.')
    parser.add_argument('-l', '--locustag', metavar='STR', type=str, default='',
        help='Input the locustag of the gene.')

    return parser


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
        sort_gff.py -- Sort genome annotation result files
attention:
        sort_gff.py genome.gff >genome.new.gff
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    sort_gff(args.gff, args.locustag)


if __name__ == "__main__":
    main()
