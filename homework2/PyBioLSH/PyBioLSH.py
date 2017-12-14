import argparse
import sys
from readfq import readfq
from minhash import BioMinHash
from lsh import BioLSH
import numpy as np


def kmer_set(seq, k):
    result = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        result.add(kmer.encode('utf8'))

    return result


def kmer_minhash(seq, k, num_perm, rc=False, debug=False):
    result = BioMinHash(num_perm=num_perm)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        result.update(kmer.encode('utf8'), rc=rc)
    if debug:
        result.output()
    return result


def main(debug=False):
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="path for sequence pasta file", type=str)
    parser.add_argument("k", help="size of kmer", type=int)
    parser.add_argument("perm", help="number of permutation", type=int)
    parser.add_argument("threshold", help="threshold to report similar", type=float)
    parser.add_argument("-cluster", help="generate clusters instead of find similar reads", required=False,
                        action='store_true')
    parser.add_argument("-debug", required=False, action='store_true')

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    handle = open(args.path)

    lsh = BioLSH(threshold=args.threshold, num_perm=args.perm)
    minhash_dict = {}
    minhash_rc_dict = {}

    if args.debug:
        set_dict = {}

    with lsh.insertion_session() as session:
        for name, seq, qual in readfq(handle):
            m_hash = kmer_minhash(seq, args.k, args.perm, debug=args.debug)

            """
            if name == 'r1':
                m_hash.set_hashvalues(np.array([0, 11, 22, 3, 5, 6, 9, 45, 98, 0, 1, 7]))
            elif name == 'r2':
                m_hash.set_hashvalues(np.array([11, 9, 3, 4, 98, 0, 1, 7, 23, 15, 0, 31]))
            """

            session.insert(name, m_hash)
            minhash_dict[name] = m_hash
            minhash_rc_dict[name] = kmer_minhash(seq, args.k, args.perm, rc=True)

            if args.debug:
                set_dict[name] = kmer_set(seq, args.k)

    if args.cluster:
        clusters = lsh.cluster()
        for cluster in clusters:
            print("\t".join(cluster))

    else:

        for id in minhash_dict.keys():
            fw_result = lsh.query(minhash_dict[id])
            rc_result = lsh.query(minhash_rc_dict[id])
            result = list(set(fw_result).union(set(rc_result)))
            for similar_id in result:
                if id != similar_id:
                    print("{}\t{}".format(id, similar_id))

    if args.debug:
        for id in minhash_dict.keys():
            for query_id in minhash_dict.keys():
                print()
                print("query_id: {} target_id: {}".format(query_id, id))
                result = minhash_dict[query_id].jaccard(minhash_dict[id])
                print("MinHash Estimated:", result)
                result = float(len(set_dict[query_id].intersection(set_dict[id]))) / float(
                    len(set_dict[query_id].union(set_dict[id])))
                print("Real Similarity", result)
                result = minhash_dict[query_id].hamming(minhash_dict[id])
                print("Hamming distance of signature vector:", result)


if __name__ == '__main__':
    main()
