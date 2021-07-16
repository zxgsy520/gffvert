#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import gzip
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)


__version__ = "v2.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_csv(file, sep='\t'):

    if file.endswith(".gz"):
        fg = gzip.open(file)
    else:
        fg = open(file)

    for line in fg:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)
    fg.close()


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


def to_string(attributes):

    attr = []

    for key, value in attributes.items():
        if key in "ID":
            attr.insert(0, '%s=%s' % (key, value))
        elif key in "Name":
            attr.insert(1, '%s=%s' % (key, value))
        elif key in "Parent":
            attr.insert(2, '%s=%s' % (key, value))
        else:
            attr.append('%s=%s' % (key, value))

    return ';'.join(attr)


def get_geneid(attributes):

    gid = ""
    attr = split_attr(attributes)

    if "ID" in attr:
        gid = attr["ID"]
    elif "locus_tag" in attr:
        gid = attr["locus_tag"]
    elif "Parent" in attr:
        gid = attr["Parent"]
    else:
        raise Exception("Gene has no name.")

    return gid


def get_cdsid(attributes):

    cdsid = ""
    attr = split_attr(attributes)

    if "Parent" in attr:
        cdsid = attr["Parent"]
    elif "locus_tag" in attr:
        cdsid = attr["locus_tag"]
    elif "ID" in attr:
        cdsid = attr["ID"]
    else:
        raise Exception('CDS has no name.')

    return cdsid


def read_sorted_gff(file):

    r = []
    stru_type = ["gene", "mRNA", "CDS"]

    for line in read_csv(file):
        if line[2] not in stru_type:
            continue
        r.append(line)

    r = sorted(r, key = lambda x:(x[0], int(x[3]), int(x[4])))

    return r


def cut_gff(file):

    gene_dict = collections.OrderedDict()
    mrna_dict = {}
    mrna2cdsid = {}
    cds_dict = {}
    gid = ""

    for line in read_sorted_gff(file):
        if line[2]=="gene":
            gid = get_geneid(line[-1])
            gene_dict[gid] = line
            continue
        if line[2]=="mRNA":
            gid = get_cdsid(line[-1])
            mrna_dict[gid] = line
            attr = split_attr(line[-1])
            mrna2cdsid[gid] = attr["ID"]
            continue

        gid = get_cdsid(line[-1])
        if gid not in cds_dict:
            cds_dict[gid] = []
        cds_dict[gid].append(line)

    return gene_dict, mrna_dict, mrna2cdsid, cds_dict


def sort_gff(file, locustag):

    n = 0
    gene_dict, mrna_dict, mrna2cdsid, cds_dict = cut_gff(file)

    print("##gff-version 3")
    for gid in gene_dict:
        n += 1
        line = gene_dict[gid]
        if gid not in mrna_dict:
            continue
        if locustag!='':
            ngid = '%s_%05d' % (locustag, n)
            LOG.info("%s\t%s" % (gid, ngid))
        else:
            ngid = gid
        attr = split_attr(line[-1])
        attr["ID"] = ngid
        attr["Name"] = ngid
        line[-1] = to_string(attr)

        print('\t'.join(line))
        if gid in mrna_dict:
            line = mrna_dict[gid]
        else:
            line[2] = 'mRNA'
        attr = split_attr(line[-1])
        attr["Parent"] = ngid
        attr["Name"] = ngid
        ngid = '%s.mrna1' % ngid
        attr["ID"] = ngid
        line[-1] = to_string(attr)
        print('\t'.join(line))

        cdsid = mrna2cdsid[gid]
        if cdsid not in cds_dict:
            LOG.info("Gene %s did not find CDS" % gid)
            continue
        k = 0
        for line in cds_dict[cdsid]:
            k += 1
            ncdsid = '%s.CDS%d' % (ngid, k)
            attr = split_attr(line[-1])
            attr["Parent"] = ngid
            attr["Name"] = ngid

            if len(cds_dict[cdsid])>=1:
                line[2] = "exon"
                exonid = '%s.exon%d' % (ngid, k)
                attr["ID"] = exonid
                line[-1] = to_string(attr)
                print('\t'.join(line))
            line[2] = "CDS"
            attr["ID"] = ncdsid
            line[-1] = to_string(attr)
            print('\t'.join(line))

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
        sort_gff.py genome.gff >genome.new.gff 2>gene.name.tsv
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    sort_gff(args.gff, args.locustag)


if __name__ == "__main__":
    main()
