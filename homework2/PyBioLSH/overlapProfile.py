import argparse
import sys
import os
from readfq import readfq
from minhash import BioMinHash
import collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

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
    parser.add_argument("groundtruth", help="path for sequence pasta file", type=str)
    parser.add_argument("k", help="size of kmer", type=int)
    parser.add_argument("perm", help="number of permutation", type=int)

    try:
        args = parser.parse_args()

    except:
        parser.print_help()
        sys.exit(1)

    handle = open(args.path)


    truth = collections.defaultdict(dict)
    with open(args.groundtruth) as f1:
        for line in f1:
            sp = line.split()
            query_id = sp[0]
            target_id = sp[1]
            overlap_size = int(sp[2])
            truth[query_id][target_id] = overlap_size

    minhash_dict = {}
    minhash_rc_dict = {}

    for name, seq, qual in readfq(handle):
        m_hash = kmer_minhash(seq, args.k, args.perm)
        minhash_dict[name] = m_hash
        minhash_rc_dict[name] = kmer_minhash(seq, args.k, args.perm, rc=True)

    identity_true = []
    identity_false = []
    identities = []
    jaccard_true = []
    jaccard_false = []
    jaccards = []
    for id in minhash_dict.keys():
        for query_id in minhash_dict.keys():
            if id != query_id:
                jaccard = minhash_dict[query_id].jaccard(minhash_dict[id])
                identity = minhash_dict[query_id].identity(minhash_dict[id])
                jaccards.append(jaccard)
                identities.append(identity)
                if id in truth[query_id]:
                    jaccard_true.append(jaccard)
                    identity_true.append(identity)
                else:
                    jaccard_false.append(jaccard)
                    identity_false.append(identity)

    cur_dir = os.path.dirname(args.path)
    cur_file = os.path.basename(args.path).split(".")[0]
    jaccard_file = os.path.join(cur_dir,cur_file+'_jaccard_k_{}_perm_{}.png'.format(args.k, args.perm))
    identity_file = os.path.join(cur_dir,cur_file+'_identity_k_{}_perm_{}.png'.format(args.k, args.perm))

    plt.figure()
    hist, bins = np.histogram(identities, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    hist1, bins = np.histogram(identity_true, bins=bins, normed = True)
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist1, align='center', width=width, label="True", alpha=0.7)
    hist1, bins = np.histogram(identity_false, bins=bins, normed = True)
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist1, align='center', width=width, label="False", alpha=0.7)
    plt.legend()
    plt.savefig(identity_file)

    plt.figure()
    hist, bins = np.histogram(jaccards, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    hist1, bins = np.histogram(jaccard_true, bins=bins, normed = True)
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist1, align='center', width=width, label="True", alpha=0.7)
    hist1, bins = np.histogram(jaccard_false, bins=bins, normed = True)
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist1, align='center', width=width, label="False", alpha=0.7)
    plt.legend()
    plt.savefig(jaccard_file)

if __name__ == '__main__':
    main()