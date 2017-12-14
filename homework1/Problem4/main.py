import argparse
import sys
from Bio import SeqIO
from global_align import linear_align as align
from global_align import linear_reverse as reverse


def align2fasta():
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="query fasta", type=str)
    parser.add_argument("target", help="target fasta", type=str)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    query_file = args.query
    target_file = args.target

    query_records = list(SeqIO.parse(query_file, "fasta"))
    target_records = list(SeqIO.parse(target_file, "fasta"))

    for query_rec in query_records:
        for target_rec in target_records:
            query_id = query_rec.id
            query_seq = query_rec.seq

            target_id = target_rec.id
            target_seq = target_rec.seq

            score = align(query_seq, target_seq)
            print("Compare {} with {}: align score is {}".format(query_id, target_id, str(score)))

            score = reverse(query_seq, target_seq)
            print("Compare {} with {}: reverse align score is {}".format(query_id, target_id, str(score)))


if __name__ == '__main__':
    align2fasta()
