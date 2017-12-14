import argparse
import sys
from readfq import readfq
from FMindex import FMindex

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="path for query sequence pasta file", type=str)
    parser.add_argument("target", help="path for target sequence pasta file", type=str)
    parser.add_argument("threshold", help="threshold to report overlap", type=int)
    parser.add_argument("-align", required=False, action='store_true')

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    # build index for all target sequence
    handle = open(args.target)
    target_dict = {}
    seqs = []
    names = []
    for name, seq, qual in readfq(handle):
        seqs.append(seq)
        names.append(name)
        target_dict[name] = seq
    print("build index")
    fm_index = FMindex(seqs, names)
    print("finish build")

    # test overlap
    handle = open(args.query)
    for name, seq, qual in readfq(handle):
        # print(name)

        outputs = fm_index.find_overlaps(seq, name, args.threshold)
        for output in outputs:
            if (output[0]!=output[3]):
                output.append("fw")
                print("\t".join(output))
                if args.align:
                    print_align(output, seq, target_dict[output[3]])

        outputs = fm_index.find_overlaps(revcomp(seq), name, args.threshold)
        for output in outputs:
            if (output[0] != output[3]):
                output.append("rc")
                print("\t".join(output))
                if args.align:
                    print_align(output, seq, target_dict[output[3]])

def print_align(output, query_seq, target_seq):
    query_len = len(query_seq)
    target_len = len(target_seq)
    query_start = int(output[1])
    query_end = int(output[2])
    target_start = int(output[4])
    target_end = int(output[5])
    print(query_seq)
    print(" "*(query_start+1) + target_seq)
    print()

if __name__ == '__main__':
    main()



