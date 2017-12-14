import argparse
import sys
from readfq import readfq
from FMindex import FMindex
import pickle

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="path for query sequence pasta file", type=str)
    parser.add_argument("index", help="path for bwt index", type=str)
    parser.add_argument("threshold", help="threshold to report overlap", type=int)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    # build index for all target sequence
    fw_index, rc_index = pickle.load(open(args.index, "rb"))

    # test overlap
    handle = open(args.query)
    for name, seq, qual in readfq(handle):
        outputs = fw_index.find_overlaps(seq, name, args.threshold)
        for output in outputs:
            if (output[0]!=output[3]):
                output.append("fw")
                output.append("suffix")
                print("\t".join(output))

        outputs = fw_index.find_overlaps(revcomp(seq), name, args.threshold)
        for output in outputs:
            if (output[0] != output[3]):
                output.append("rc")
                output.append("suffix")
                print("\t".join(output))

        outputs = rc_index.find_overlaps(revcomp(seq), name, args.threshold)
        for output in outputs:
            if (output[0] != output[3]):
                output.append("fw")
                output.append("prefix")
                print("\t".join(output))

        outputs = rc_index.find_overlaps(seq, name, args.threshold)
        for output in outputs:
            if (output[0] != output[3]):
                output.append("rc")
                output.append("prefix")
                print("\t".join(output))


if __name__ == '__main__':
    main()