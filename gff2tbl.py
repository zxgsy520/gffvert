#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO) #初始化时，如果没指定level，那么level的默认级别为WARNING

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang,Shuying Deng",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    """读取表格文件"""
    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line or line.startswith("#"): #跳过空行和首字母为#的行
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


def split_attr(attributes, sep="="):

    """分割gff文件中的最后一列"""
    r = collections.OrderedDict()
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if sep not in content:
            LOG.info("%r is not a good formated attribute: no tag!" % content)
            continue
        tag, value = content.strip().split(sep, 1)
        r[tag] = value.strip('"')

    return r


def read_gff2gene(file):

    gene = []
    mrna_dict = {}
    cds_dict = {}
    geneid = ""
    parent = ""

    for line in read_tsv(file, "\t"):
        if line[2] == "gene":
            gene.append(line)
            attr = split_attr(line[-1])
            geneid = attr["ID"]
            continue
        if line[2] == "mRNA":
            attr = split_attr(line[-1])
            if geneid != attr["Parent"]:
                LOG.info("Gene %s has no mRNA information" % geneid)
                break
            parent = attr["ID"]
            if geneid not in mrna_dict:
                mrna_dict[geneid] = line

        if line[2] == "CDS":
            attr = split_attr(line[-1])
            if parent != attr["Parent"]:
                LOG.info("Gene %s has no CDS information" % geneid)
                break
            if geneid not in cds_dict:
                cds_dict[geneid] = []
            cds_dict[geneid].append(line)

    return gene, mrna_dict, cds_dict


def gff2tbl(gff, function):

    r = read_function(function)
    gene, mrna_dict, cds_dict = read_gff2gene(gff)
    seqid = ""

    for line in sorted(gene, key = lambda x:(x[0], int(x[3]), int(x[4]))):
        attr = split_attr(line[-1])
        geneid = attr["ID"]
        if geneid not in cds_dict:
            LOG.info("Gene %s has no CDS information" % geneid)
            break
        if geneid not in mrna_dict:
           LOG.info("Gene %s has no mRNA information" % geneid)
           break
        if line[0] != seqid:
            print(">Feature %s" % line[0])
            seqid = line[0]
        strand = line[6]
        mrna = mrna_dict[geneid]
        mrnaid = split_attr(mrna[-1])["ID"]

        gene_nane = ""
        func = "hypothetical protein"
        if mrnaid in r:
            gene_nane = r[mrnaid][0]
            func = r[mrnaid][1]
        if gene_nane:
            gene_nane = "\n\t\tgene\t%s" % gene_nane

        print("{start}\t{end}\tgene{gene}\n\t\tlocus_tag\t{locus_tag}".format(start=line[3],
            end=line[4], gene=gene_nane, locus_tag=geneid))
        temp = []
        for line in cds_dict[geneid]:
            if line[6] == "-":
                temp.append((line[4], line[3]))
            else:
                temp.append((line[3], line[4]))
        #if strand == "-":
        #    temp = temp[::-1]
        print("%s\t%s\tCDS" % (temp[0][0], temp[0][1]))
        for i in temp[1::]:
            print("%s\t%s" % (i[0], i[1]))
        print("\t\tproduct\t%s\n\t\tprotein_id\t%s" % (func, mrnaid))
        

    return 0


def add_help_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input gff file")
    parser.add_argument("-f", "--function", metavar="FILE", type=str, required=True,
        help="Input gene function annotation results, merge.function.xls.")

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
    gff2tbl.py: Convert gff files to table format files
For exmple:
    gff2tbl.py prefix.gff3 --function merge.function.xls >prefix.tbl 

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    gff2tbl(args.input, args.function)


if __name__ == "__main__":

    main()

