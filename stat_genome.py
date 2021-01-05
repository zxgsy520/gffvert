#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "v1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "1131978210@qq.com"
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
        if type(line) == type(b''):
            line = line.decode('utf-8')

        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            if seq:
                yield seq.split('\n')
            seq = "%s\n" % line
        else:
            seq += line
    if seq:
        yield seq.split('\n')
    fp.close()


def read_fastq(file):
    '''Read fastq file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq.append(line.strip("@"))
            continue
        seq.append(line)

        if len(seq)==4:
            yield seq
            seq = []
    fp.close()


def scaff2contig(seq):

    seq = seq.upper()
    r = []
    gaplen = seq.count('N')

    for i in seq.split("N"):
        if len(i) == 0:
            continue
        r.append(i)

    return r, len(r)-1, gaplen


def read_seq(file):

    scaff = {}
    contig = {}

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        seq = line[1].upper()
        seqs, gap, gaplen = scaff2contig(seq)
        scaff[line[0]] = [seq, gap, gaplen]

        if len(seqs)==1:
            contig[line[0]] = [seqs[0], 0, 0]
        else:
            n = 0
            for i in seqs:
                n += 1
                contig["%s_%s" % (line[0], n)] = [i, 0, 0]

    return scaff, contig


def stat_length(data):

    r = {}

    for i in data:
        line = data[i]
        seq, gap, gaplen = line
        seqlen = len(seq)

        if seqlen >= 1000:
            if "Length>=1kb" not in r:
                r["Length>=1kb"] = [0, 0, 0, 0]
            r["Length>=1kb"][0] += len(seq)
            r["Length>=1kb"][1] += 1
            r["Length>=1kb"][2] += gaplen
            r["Length>=1kb"][3] += gap
        if seqlen >= 2000:
            if "Length>=2kb" not in r:
                r["Length>=2kb"] = [0, 0, 0, 0]
            r["Length>=2kb"][0] += len(seq)
            r["Length>=2kb"][1] += 1
            r["Length>=2kb"][2] += gaplen
            r["Length>=2kb"][3] += gap
        if seqlen >= 5000:
            if "Length>=5kb" not in r:
                r["Length>=5kb"] = [0, 0, 0, 0]
            r["Length>=5kb"][0] += len(seq)
            r["Length>=5kb"][1] += 1
            r["Length>=5kb"][2] += gaplen
            r["Length>=5kb"][3] += gap
        if "Total" not in r:
            r["Total"] = [0, 0, 0, 0]
        r["Total"][0] += len(seq)
        r["Total"][1] += 1
        r["Total"][2] += gaplen
        r["Total"][3] += gap

    return r


def nx0(sumlen, data, x="N50"):

    dlen = {"Longest":0, "N50":0.5, "N60":0.6, "N70":0.7, "N80":0.8, "N90": 0.9}
    r = {}
    contlen = 0

    for line in sorted(data.items(), key=lambda x:x[1][0], reverse=True):
        seqlen, gap, gaplen = line[1]
        contlen += seqlen
        if x not in r:
            r[x] = [0, 0, 0, 0]
        r[x][1] += 1
        r[x][2] += gaplen
        r[x][3] += gap

        if contlen >= sumlen*dlen[x]:
            if seqlen >= r[x][0]:
                r[x][0] = seqlen
            else:
                r[x][1] += -1
                r[x][2] += -gaplen
                r[x][3] += -gap
                break
    return r


def stat_quality(data):

    sumlen = 0
    for i in data:
        sumlen += len(data[i][0])
        data[i][0] = len(data[i][0])

    r = {}
    dlen = ["Longest", "N50", "N60", "N70", "N80", "N90"]

    for x in dlen:
        r.update(nx0(sumlen, data, x))

    return r


def format_dict(data):

    r = {}

    for i in data:
        line = []
        for j in data[i]:
            try:
                j = "{0:,}".format(int(j))
            except:
                j = "{0:,}".format(float(j))
            else:
                j = j
            line.append(j)
        r[i] = line

    return r


def stat_genome(file):

    scaff, contig = read_seq(file)
    slen = format_dict(stat_length(scaff))
    snx = format_dict(stat_quality(scaff))

    if scaff != contig:
        clen = format_dict(stat_length(contig))
        cnx = format_dict(stat_quality(contig))
    else:
        clen = slen
        cnx = snx

    print("""\
#StatType\tContigLength\tContigNumber\tScaffoldLength\tScaffoldNumber\tGapLength\tGapNumber
N50\t{0}\t{1}
N60\t{2}\t{3}
N70\t{4}\t{5}
N80\t{6}\t{7}
N90\t{8}\t{9}
Longest\t{10}\t{11}
Total\t{12}\t{13}
Length>=1kb\t{14}\t{15}
Length>=2kb\t{16}\t{17}
Length>=5kb\t{18}\t{19}\
""".format("\t".join(cnx["N50"][0:2]), "\t".join(snx["N50"]),
    "\t".join(cnx["N60"][0:2]), "\t".join(snx["N60"]),
    "\t".join(cnx["N70"][0:2]), "\t".join(snx["N70"]),
    "\t".join(cnx["N80"][0:2]), "\t".join(snx["N80"]),
    "\t".join(cnx["N90"][0:2]), "\t".join(snx["N90"]),
    "\t".join(cnx["Longest"][0:2]), "\t".join(snx["Longest"]),
    "\t".join(clen["Total"][0:2]), "\t".join(slen["Total"]),
    "\t".join(clen["Length>=1kb"][0:2]), "\t".join(slen["Length>=1kb"]),
    "\t".join(clen["Length>=2kb"][0:2]), "\t".join(slen["Length>=2kb"]),
    "\t".join(clen["Length>=5kb"][0:2]), "\t".join(slen["Length>=5kb"])
    ))


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input the genome file.')

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
    stat_genome -- Statistical genome length and base contentStat N50-N90 value of contigs and scaffolds.

attention:
    stat_genome --fasta genome.fasta
''')
    args = add_hlep_args(parser).parse_args()
    stat_genome(args.input)


if __name__ == "__main__":
    main()
