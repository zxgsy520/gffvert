#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import logging
import argparse


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


def sorted_pos(start, end):

    if start >end:
        start, end = end, start

    return start, end


def convert_pos(rbeg, qbeg, qend, pos, direct='+'):

    npos = 0

    if direct=='+':
        npos = (pos-rbeg)+qbeg
    elif direct=='-':
        npos = qend-(pos-rbeg)
    else:
        raise Exception('%s is not in [-, +]' % direct)
    if  npos <= 0:
        raise Exception('Please check:rbeg={rbeg}, qbeg={qbeg}, qend={qend}, pos={pos}, direct={direct}'.format(**locals()))

    return npos


def read_bed(file):

    data = {}

    for line in read_csv(file):
        if line[0] not in data:
            data[line[0]] = []
        rstart, rend = sorted_pos(int(line[1]), int(line[2]))
        qstart, qend = sorted_pos(int(line[6]), int(line[7]))
        data[line[0]].append([[rstart, rend], [qstart, qend], line[4], line[5]])

    return data


def judge_att(dlist, rstart, rend):

    nlist = False

    for line in dlist:
        rpot = line[0]
        if rstart>=rpot[0] and rend<=rpot[1]:
            nlist = line
            break

    return nlist


def convert_position_gff(gff, bed):

    data = read_bed(bed)

    for line in read_csv(gff):
        if line[0] not in data:
            LOG.debug('Did not change the %s coordinate of the sequence' % line[0])
            print('\t'.join(line))
            continue
        rstart, rend = sorted_pos(int(line[3]), int(line[4]))
        nlist = judge_att(data[line[0]], rstart, rend)

        if not nlist:
            continue
        rbeg = nlist[0][0]
        qbeg, qend = nlist[1]
        direct = nlist[-1]
        nstart = convert_pos(rbeg, qbeg, qend, rstart, direct)
        nend = convert_pos(rbeg, qbeg, qend, rend, direct)
        ndirect = line[6]

        if direct=='-':
            if ndirect=='+':
                ndirect = '-'
            else:
                ndirect = '+'

        nstart, nend = sorted_pos(nstart, nend)

        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
            nlist[-2], line[1], line[2], nstart, nend, line[5], ndirect, line[7], line[8]))


def add_help_args(parser):

    parser.add_argument('gff',
        help='Input genome annotation gff file.')
    parser.add_argument('-b', '--bed', metavar='STR', type=str, default='protein',
        help='Input the bed file mounted by Hi-C.')

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
        convert_position_gff.py -- Transform gene coordinates
attention:
        convert_position_gff.py genome.gff --bed hic.bed  >genome.new.gff
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    convert_position_gff(args.gff, args.bed)


if __name__ == "__main__":
    main()
