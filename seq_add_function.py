#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.2"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep="\t"):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fh.close()


def get_best_gene(gstr):
 
    r = ""
    for i in gstr.split(","):
        if not i:
            continue
        if len(i) <= 4 or len(i) >5:
            continue
        r = i

    return r


def read_function(file):

    r = {}

    for line in read_tsv(file, sep="\t"):
        gene_id = line[0]
        gene = line[2]
        if gene == "-":
            gene = ""
        gene = gene.replace(" ", "").replace("-", "").replace(",,", ",").strip(",")
        gene = get_best_gene(gene)
        func = line[3]
        if func == "-":
            func = line[6]
        if func == "-":
            func = line[10]
        if func == "-":
            func = line[13]
        if func == "-":
            func = line[16]
        if func == "-":
            func = "hypothetical protein"
        r[gene_id] = [gene, func]

    return r


def seq_add_function(file, function, origin="Unknown"):

    r = read_function(function)

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".faa") or file.endswith(".fna"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()
        if line.startswith(">"):
            seqid = line.strip(">").split()[0]
            temp = ""
            if seqid in r:
                gene, func = r[seqid]
                if gene != "":
                    temp = "[gene=%s]" % gene
                if func != "":
                    temp = "%s [protein=%s]" % (temp, func)
                if origin:
                    temp = "%s [organism=%s]" % (temp, origin)
                temp = temp.replace("  ", "").strip()
            line = "%s %s" % (line.split()[0], temp)
            line = line.strip()
                  
        print(line)

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input sequence")
    parser.add_argument("-f", "--function", metavar="FILE", type=str, required=True,
        help="Input gene function annotation results, PVgenome.merge.function.xls.")
    parser.add_argument("-or", "--origin", metavar="FILE", type=str, default="Unknown",
        help="Input origin species information")

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
    seq_add_function.py: Add information to the sequence
For exmple:
    seq_add_function.py protein.fasta --function merge.function.xls --origin "Fonsecaea pedrosoi" > prefix.pep

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    seq_add_function(args.input, args.function, args.origin)


if __name__ == "__main__":

    main()
