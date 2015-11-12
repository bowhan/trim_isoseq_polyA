#!/usr/bin/env python
"""
this script filter the fasta sequences by their percentage of A
usage:
    -i input.fasta -t A -p 0.75 > reads_of_insert.rna.fasta
"""


def filter(filename, target, percentage):
    i = 0
    j = 0
    from Bio import SeqIO
    percentage = float(percentage)
    for fa in SeqIO.parse(open(filename), "fasta"):
        if (1.0 * fa.seq.count(target)) / len(fa.seq) > percentage:
            print ">%s\n%s" % (fa.id, fa.seq)
            i += 1
        j += 1
    return (i, j)


if __name__ == "__main__":
    import argparse, sys

    parser = argparse.ArgumentParser(description='this script filter the fasta sequences by their percentage of A')
    parser.add_argument('-i', required=True, action="store", type=str, dest="input",
                        help="fasta file storing the reads")
    parser.add_argument('-t', required=False, action="store", type=str, dest="target", default='A',
                        help="the target nucleotide")
    parser.add_argument('-p', required=False, action="store", type=float, dest="threshold", default=0.75,
                        help="the percentage threshold; percentage lower than this will NOT be printed")
    args = parser.parse_args()
    (trimmed, total) = filter(args.input, args.target, args.threshold)
    print >> sys.stderr, "%d\t%d" % (total, trimmed)
