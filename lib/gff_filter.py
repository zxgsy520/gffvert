#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "v1.1.3"
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


def read_gff(file):

    r = []
    for line in read_tsv(file, "\t"):
        if line[2] in ["mRNA", "transcript"]:
            if len(r) >= 1:
                yield r
            r = []
        r.append(line)
    if len(r) >= 1:
        yield r


def filter_gff(data, cds=3, min_cds=900, max_cds=2000, cds_exon=0.6, intron=80):

    elen = 0
    cdslen = 0
    cds_number = 0
    r = data

    for line in data:
        start, end = int(line[3]), int(line[4])
        if start >= end:
            start, end = end, start
        if line[2] in ["mRNA", "transcript"]:
            instart = start
            continue
        if line[2] == "exon":
            if start-instart > intron:
                LOG.info('%s\t%s\t%s' % (line[-1], start, instart))
                r = []
                break
            instart = end
            elen += (end-start+1)
            continue
        if line[2] == "CDS":
            cds_number += 1
            cdslen += (end-start+1)
    if cds_number>=3 and (cdslen>=min_cds and cdslen<=max_cds) and cdslen*1.0/elen>=cds_exon:
        pass
    else:
        LOG.info('%s\t%s\t%s\t%s' % (line[-1], cds_number, cdslen, elen))
        r = []

    return r


def gff_filter(file, cds=3, min_cds=900, max_cds=2000, cds_exon=0.6, intron=80):

    for r in read_gff(file):
        r = filter_gff(r, cds, min_cds, max_cds, cds_exon, intron)
        if len(r) <=1:
            continue
        for line in r:
            print('\t'.join(line))
    return 0


def add_help_args(parser):

    parser.add_argument('gff',
        help='Input genome annotation gff file.')
    parser.add_argument('-c', '--cds', metavar='INT', type=int, default=3,
        help='Set the minimum number of CDS, default=3.')
    parser.add_argument('--min_cds', metavar='INT', type=int, default=900,
        help='Set the minimum CDS length (bp), default=900.')
    parser.add_argument('--max_cds', metavar='INT', type=int, default=5000,
        help='Set the maximum CDS length (bp), default=5000.')
    parser.add_argument('--cds_exon', metavar='FLOAT', type=float, default=0.6,
        help='Set the minimum ratio of CDS to exon, default=0.6.')
    parser.add_argument('--intron', metavar='INT', type=int, default=100,
        help='Set the maximum length of intron, default=80.')

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
        gff_filter.py -- Filter the results of PASA for training.
attention:
        gff_filter.py trainset.gtf >trainset_new.gft
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    gff_filter(args.gff, args.cds, args.min_cds, args.max_cds, args.cds_exon, args.intron)


if __name__ == "__main__":
    main()
