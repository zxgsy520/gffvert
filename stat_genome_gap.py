#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "v1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip("@").split()[0]
            continue
        if line.startswith('@'):
            line = line.strip("@").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq = "%s\n" % line
        else:
            seq += "%s\n" % line

    if len(seq.split('\n'))==5:
        seq = seq.split('\n')
        yield seq[0], seq[1]


def stat_gap(seq):

    n = 0
    start = 0
    tep = ''
    site = 0
    data = {}

    for i in seq:
        site +=1

        if i!="N":
            if tep=="" or tep!="N":
                tep = i
                continue
            else:
                tep = i
                n += 1
                data[n] = [start, site-1]
        else:
            if tep=="" or tep!="N":
                tep = i
                start = site
            else:
                tep = i
                continue
    if i=="N":
        n += 1
        data[n] = [start, site-1]

    return data


def stat_genome_gap(file):

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)
    
    print('#Seq id\tGap number\tGap start\tGap end\Gap length')
    for seqid, seq in fh:
        seq = seq.upper()
        if "N" not in seq:
            continue
        gap_dict = stat_gap(seq)

        for i in gap_dict:
            gaplen = gap_dict[i][1]-gap_dict[i][0]+1
            print('{0}\t{1}\t{2:,}\t{3:,}\t{4:,}'.format(seqid, i, gap_dict[i][0], gap_dict[i][1], gaplen))


def add_help_args(parser):

    parser.add_argument('genome',
        help='Input genome file.')

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
        stat_genome_gap -- Count the number of gaps in the assembled genome.
attention:
        stat_genome_gap genome.fasta >stat.gap.txt
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    stat_genome_gap(args.genome)


if __name__ == "__main__":
    main()
