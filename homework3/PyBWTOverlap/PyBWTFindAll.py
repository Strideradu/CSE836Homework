"""
Given an input read, these program will keep find all possible alignment
"""

import argparse
import sys
from readfq import readfq
from FMindex import FMindex
import pickle
from Bio import SeqIO

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("seed", help="path for query sequence pasta file", type=str)
    parser.add_argument("data", help="path for all reads", type=str)
    parser.add_argument("index", help="path for bwt index", type=str)
    parser.add_argument("threshold", help="threshold to report overlap", type=int)
    parser.add_argument("output", help="path for output fasta", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    # build index for all target sequence
    seq_dict = SeqIO.index(args.data, 'fasta')
    sets = []
    seeds = set()
    total = set()
    handle = open(args.seed)
    for name, seq, qual in readfq(handle):
        seeds.add(name)
    sets.append(seeds)
    total = total.union(seeds)

    fw_index, rc_index = pickle.load(open(args.index, "rb"))

    new_set = seeds
    while len(new_set) > 0:
        prev_set = new_set
        new_set = set()
        for name in prev_set:

            seq = str(seq_dict[name].seq)
            outputs = fw_index.find_overlaps(seq, name, args.threshold)
            for output in outputs:
                if output[3] not in total:
                    new_set.add(output[3])

            outputs = fw_index.find_overlaps(revcomp(seq), name, args.threshold)
            for output in outputs:
                if output[3] not in total:
                    new_set.add(output[3])


            outputs = rc_index.find_overlaps(revcomp(seq), name, args.threshold)
            for output in outputs:
                if output[3] not in total:
                    new_set.add(output[3])

            outputs = rc_index.find_overlaps(seq, name, args.threshold)
            for output in outputs:
                if output[3] not in total:
                    new_set.add(output[3])


        sets.append(new_set)
        total = total.union(new_set)
        print("Find one set of reads {}".format(len(new_set)))
        print('Total number of reads is {}'.format(len(total)))

    outputs = []
    for seq_id in total:
        outputs.append(seq_dict[seq_id])

    SeqIO.write(outputs, open(args.output, 'w'), 'fasta')



if __name__ == '__main__':
    main()