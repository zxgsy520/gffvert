#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "v1.1.4"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def split_attr(attributes):

    r = OrderedDict()
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if "=" not in content:
            LOG.info("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def to_string(data):

    attr = []
    for key, value in data.items():
        if key in "ID":
            attr.insert(0, '%s=%s' % (key, value))
        else:
            attr.append('%s=%s' % (key, value))

    return ';'.join(attr)


def get_mrna2gene(exon):

    exon[2] = "gene"
    exon[-1]["ID"] = exon[-1]["Parent"]

    return exon


def read_gff(file):

    data = OrderedDict()
    mrna2gene = {}

    for line in read_tsv(file, "\t"):
        line[-1] = split_attr(line[-1])

        if line[2] == "gene":
            geneid = line[-1]["ID"]
            data[geneid] = {"mRNA": {"CDS": [], "exon": []}, "gene": line}
            continue

        if line[2] == "mRNA":
            geneid = line[-1]["Parent"]
            mrnaid = line[-1]["ID"]
            mrna2gene[mrnaid] = geneid
            if geneid not in data:
                gene = get_mrna2gene(line)
                geneid = gene[-1]["ID"]
                data[geneid] = {"mRNA": {"CDS": [], "exon": []}, "gene": line}
            data[geneid]["mRNA"]["feature"] = line
            continue

        if line[2] in ["exon", "CDS"]:
            mrnaid = line[-1]["Parent"]
            if mrnaid not in mrna2gene:
                continue
            geneid = mrna2gene[mrnaid]
            data[geneid]["mRNA"][line[2]].append(line)

    return data


def rename_evm_gff(file, gene_format="CgT%sg%05d0"):

    data = read_gff(file)
    n = 0
    chr = ""
    k = 0

    for geneid in data:
        line = data[geneid]["gene"]
        if line[0] != chr:
            chr = line[0]
            k = 0
            n += 1
        k += 1
        new_geneid = gene_format % (n, k)
        line[-1]["ID"] = new_geneid
        line[-1]["Name"] = new_geneid
        line[-1] = to_string(line[-1])
        print("\t".join(line))

        line = data[geneid]["mRNA"]["feature"]
        new_mrnaid = "%s.1" % new_geneid
        line[-1]["ID"] = new_mrnaid
        line[-1]["Name"] = new_mrnaid
        line[-1]["Parent"] = new_geneid
        line[-1] = to_string(line[-1])
        print("\t".join(line))

        cds = data[geneid]["mRNA"]["CDS"]
        j = 0
        for line in data[geneid]["mRNA"]["exon"]:
            j += 1
            exonid = "%s.exon%s" % (new_mrnaid, j)
            line[-1]["ID"] = exonid
            line[-1]["Parent"] = new_mrnaid
            line[-1] = to_string(line[-1])
            print("\t".join(line))

            line = cds[j-1]
            cdsid = "%s.cds%s" % (new_mrnaid, j)
            line[-1]["ID"] = cdsid
            line[-1]["Parent"] = new_mrnaid
            line[-1] = to_string(line[-1])
            print("\t".join(line))

    return 0


def add_help_args(parser):

    parser.add_argument("gff",  metavar="FILE", type=str,
        help="Input genome annotation gff file.")
    parser.add_argument("-gf", "--gene_format", metavar="STR", type=str,
        help="Input the format for renaming genes.")

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
        rename_evm_gff.py -- Duplicate the name of the gene in the gff file.
attention:
        rename_evm_gff.py genome.gff  >genome.new.gff
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    rename_evm_gff(args.gff, args.gene_format)


if __name__ == "__main__":
    main()
