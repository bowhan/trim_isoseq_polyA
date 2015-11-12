#!/usr/bin/env python
"""
this script extracts the insert + polyA sequence from "reads_of_insert.fasta"
with 5' and 3' end information from "isoseq_draft.primer_info.csv"
usage:
    -c isoseq_draft.primer_info.csv -i reads_of_insert.fasta > reads_of_insert.rna.fasta
"""


def parse_primer_info(primer_info):
    from collections import defaultdict, namedtuple
    pos_line = namedtuple("primer_info_tokens",
                          "id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera")
    pos_map = defaultdict(pos_line)
    for line in open(primer_info):
        tokens = line.split(",")
        name = '/'.join(tokens[0].split('/')[0:2]) + "/ccs"
        pos_map[name] = pos_line._make(tokens)
    return pos_map


def extract_fa_seq(fa_file, pos_map):
    from Bio import SeqIO
    for fa in SeqIO.parse(fa_file, "fasta"):
        pos_info = pos_map[fa.id] if fa.id in pos_map else None
        if pos_info \
                and pos_info.polyAseen == '1'  \
                and pos_info.threeseen == '1'  \
                and pos_info.fiveend   != 'NA' \
                and pos_info.threeend  != 'NA' :
            if pos_info.strand == '+':
                print ">%s\n%s" % (fa.id, fa.seq[int(pos_info.fiveend) - 1: int(
                    pos_info.polyAend) + 4])
                # TODO: figure out why there is offset of 4
            else:
                print ">%s\n%s" % (
                    fa.id, fa.seq.reverse_complement()[int(pos_info.fiveend) - 1: int(pos_info.threeend) + 4])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='this scrpit extracts polyA sequence')
    parser.add_argument('-c', required=True, action="store", type=str, dest="primer_info",
                        help="isoseq_draft.primer_info.csv which stores the coordinate information")
    parser.add_argument('-i', required=True, action="store", type=str, dest="input",
                        help="reads_of_insert.fasta which stores the reads of insert files")
    args = parser.parse_args()
    pos_map = parse_primer_info(args.primer_info)
    extract_fa_seq(args.input, pos_map)
