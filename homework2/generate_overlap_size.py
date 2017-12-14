# from alignment file from PBsim to generate overlap dict
from intervaltree import Interval, IntervalTree
from Bio import AlignIO

maf_file = "D:/Data/hmmer for pacbio/Pacbio simulate/Ecoli_Pacbio_simulate_20X.maf"
# maf_file = "C:/Users/dunan/Documents/GitHub/CSE836_HW/homework2/P4_30X_0001.maf"
tree = IntervalTree()
seq_dict = {}
for multiple_alignment in AlignIO.parse(maf_file, "maf"):
    multiple_alignment = list(multiple_alignment)
    id = multiple_alignment[1].id
    start = multiple_alignment[0].annotations["start"]
    end = start + multiple_alignment[0].annotations["size"]
    tree[start:end] = id
    seq_dict[id] = (start, end)

with open("C:/Users/dunan/Documents/GitHub/CSE836_HW/homework2/Ecoli_20X_overlap.txt", "w") as fout:
    seq_list = list(seq_dict.keys())
    overlap_dict = {}
    for seq_id in seq_list:
        #print(seq_id)
        overlap_list = list(tree.search(seq_dict[seq_id][0], seq_dict[seq_id][1]))
        for overlap_rec in overlap_list:
            if overlap_rec.data != seq_id:
                target_id = overlap_rec.data
                x = range(seq_dict[seq_id][0], seq_dict[seq_id][1])
                y = range(overlap_rec.begin, overlap_rec.end)
                overlap_len = len(set(x) & set(y))
                if overlap_len > 100:
                    overlap_dict[(seq_id, target_id)] = True
                    print("{}\t{}\t{}".format(seq_id,target_id, overlap_len), file=fout)


