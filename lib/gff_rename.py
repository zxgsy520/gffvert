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


def split_attr(attributes):

    r = OrderedDict()
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


def to_string(data):

    attr = []
    for key, value in data.items():
        if key in "ID":
            attr.insert(0, '%s=%s' % (key, value))
        else:
            attr.append('%s=%s' % (key, value))

    return ';'.join(attr)


def gff_rename(file, prefix):

    r = []
    gn = 0
    for line in read_tsv(file, "\t"):
        if line[2] == "gene":
            if len(r) >= 2:
                yield r
            r = []
            en = 0
            gn += 1
            gi = "%s_%s_%s" % (prefix, line[0], gn)
        attr = split_attr(line[-1])
        if "ID" in attr:
            attr["ID"] = gi
        if "Parent" in attr:
            attr["Parent"] = gi
        if "Name" in attr:
            attr["Name"] = "prediction_%s_%s" % (line[0], gn)
        if line[2] == "exon":
            en += 1
            ei = "%s.exon%s" % (gi, en)
            if "ID" in attr:
                attr["ID"] = ei
        line[-1] = to_string(attr)
        r.append(line)
    if len(r) >= 2:
        yield r


def out_gff(file, prefix):

    for r in gff_rename(file, prefix):
        for line in r:
            print('\t'.join(line))
        print("")

def add_help_args(parser):

    parser.add_argument('gff',
        help='Input genome annotation gff file.')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str, default='EVM',
        help='Input the prefix of the gene name, default=EVM.')

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
        gff_rename.py -- Duplicate the name of the gene in the gff file.
attention:
        gff_rename.py genome.gff >genome.new.gff
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    out_gff(args.gff, args.prefix)


if __name__ == "__main__":
    main()
