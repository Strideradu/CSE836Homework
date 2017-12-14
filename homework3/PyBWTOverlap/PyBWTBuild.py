import argparse
import sys
from readfq import readfq
from FMindex import FMindex
import pickle

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="path for the fasta file to build bwt", type=str)
    parser.add_argument("save", help="path to save bwt", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    handle = open(args.input)

    seqs = []
    rc_seqs = []
    names = []
    for name, seq, qual in readfq(handle):
        seqs.append(seq)
        rc_seqs.append(revcomp(seq))
        names.append(name)
    print("build forward index")
    fm_index = FMindex(seqs, names)
    print("finish build")

    print("build reverse index")
    rc_index = FMindex(rc_seqs, names)
    print("finish build")
    index = [fm_index, rc_index]
    pickle.dump(index, open(args.save, "wb"))

if __name__ == '__main__':
    main()

