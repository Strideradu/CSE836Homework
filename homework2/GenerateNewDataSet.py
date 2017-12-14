from intervaltree import Interval, IntervalTree
from Bio import AlignIO
from Bio import SeqIO

num_seqs = 350
upper_overlap = 3000
lower_overlap = 500
read_length = 3000

longer = set()
shorter = set()

maf_file = "/mnt/home/dunan/RawReads/EColi_MG1655_Pacbio/P2C4_30X_simulated/P4_30X_0001.maf"
overlap_txt = "/mnt/home/dunan/Learn/Class/CSE836/CSE836_HW/homework2/Ecoli_30X_overlap.txt"
fasta_file = "/mnt/home/dunan/RawReads/EColi_MG1655_Pacbio/P2C4_30X_simulated/Ecoli_simulated_30X_C4.fasta"
save_fasta = "/mnt/home/dunan/Learn/Class/CSE836/CSE836_HW/homework2/Ecoli_30x_overlap_{}_{}_cutoff.fasta".format(upper_overlap, lower_overlap)
new_overlap_txt = "/mnt/home/dunan/Learn/Class/CSE836/CSE836_HW/homework2/Ecoli_30X_overlap_{}_{}_cutoff.txt".format(upper_overlap, lower_overlap)
record_dict = SeqIO.index(fasta_file, "fasta")

with open(overlap_txt) as f:
    for line in f:
        line = line.rstrip()
        if line != "" and line[0] != "#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            if len(record_dict[query_id])>=read_length and len(record_dict[target_id])>=read_length:
                overlap_size = int(line_sp[2])
                if overlap_size >= upper_overlap:
                    if len(longer) <= num_seqs:
                        longer.add(query_id)
                        longer.add(target_id)
                if overlap_size <= lower_overlap:
                    if len(shorter) <= num_seqs:
                        shorter.add(query_id)
                        shorter.add(target_id)

        if len(longer) > num_seqs and len(shorter) > num_seqs:
            break
seq_list = list(longer.union(shorter))

tree = IntervalTree()
seq_dict = {}
for multiple_alignment in AlignIO.parse(maf_file, "maf"):
    multiple_alignment = list(multiple_alignment)
    id = multiple_alignment[1].id
    start = multiple_alignment[0].annotations["start"]
    end = start + multiple_alignment[0].annotations["size"]
    tree[start:end] = id
    seq_dict[id] = (start, end)

fasta_output = []
with open(new_overlap_txt, "w") as fout:    
    overlap_dict = {}
    for seq_id in seq_list:
        fasta_output.append(record_dict[seq_id])
        #print(seq_id)
        overlap_list = list(tree.search(seq_dict[seq_id][0], seq_dict[seq_id][1]))
        for overlap_rec in overlap_list:
            if overlap_rec.data != seq_id:
                target_id = overlap_rec.data
                x = range(seq_dict[seq_id][0], seq_dict[seq_id][1])
                y = range(overlap_rec.begin, overlap_rec.end)
                overlap_len = len(set(x) & set(y))
                if overlap_len >= upper_overlap:
                    overlap_dict[(seq_id, target_id)] = True
                    print("{}\t{}\t{}".format(seq_id,target_id, overlap_len), file=fout)

SeqIO.write(fasta_output,save_fasta,"fasta")