from Bio import SeqIO
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sam", help="path for sam file", type=str)
    parser.add_argument("data", help="path for reads", type=str)
    parser.add_argument("out", help="path for output", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)
    seq_dict = SeqIO.index(args.data, 'fasta')
    output = []
    with open(args.sam) as f:
        for line in f:
            if line and line[0] != '@':
                sp = line.split()
                if sp[2]!='*':
                    output.append(seq_dict[sp[0]])

    SeqIO.write(output, open(args.out, 'w'), 'fasta')

if __name__ == '__main__':
    main()

