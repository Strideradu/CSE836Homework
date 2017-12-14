from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

fig_path = "/mnt/home/dunan/Learn/Class/CSE836/Ecoli_30x_overlap_3000_500_cutoff_k10_1024_0p04.png"
result = "/mnt/home/dunan/Learn/Class/CSE836/Ecoli_30x_overlap_3000_500_cutoff_k10_1024_0p04.out"

size = 0
for rec in SeqIO.parse("/mnt/home/dunan/Learn/Class/CSE836/CSE836_HW/homework2/Ecoli_30x_overlap_3000_500_cutoff.fasta", "fasta"):    
    size += 1

total_pairs = size*(size - 1)/2

expect = 0
found_overlap = {}
with open("/mnt/home/dunan/Learn/Class/CSE836/CSE836_HW/homework2/Ecoli_30X_overlap_3000_500_cutoff.txt") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0] != "#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            overlap_size = int(line_sp[2])
            # print(query_id, target_id)
            if found_overlap.get((query_id, target_id)) is None and found_overlap.get((target_id, query_id)) is None:
                expect += 1
                found_overlap[(query_id, target_id)] = overlap_size

all_overlap = {}
with open("/mnt/home/dunan/Learn/Class/CSE836/CSE836_HW/homework2/Ecoli_30X_overlap.txt") as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0] != "#":
            line_sp = line.split("\t")
            query_id = line_sp[0]
            target_id = line_sp[1]
            overlap_size = int(line_sp[2])
            # print(query_id, target_id)
            if all_overlap.get((query_id, target_id)) is None and all_overlap.get((target_id, query_id)) is None:
                all_overlap[(query_id, target_id)] = overlap_size

num_found = 0
true_align = 0
found = {}

find = []
notfind = []
both = []

with open(result) as f1:
    for line in f1:
        line = line.rstrip()
        if line != "" and line[0]!="o":
            line_sp = line.split("\t")
            # print(line_sp)
            query_id = line_sp[0].strip()            
            target_id = line_sp[1].strip()
            if query_id!=target_id:
                if (found.get((query_id, target_id), False) is False) and (found.get((target_id, query_id), False) is False):

                    num_found += 1
                    found[(query_id, target_id)] = True
                    if found_overlap.get((query_id, target_id), False) or found_overlap.get((target_id, query_id), False):
                        true_align += 1
                        if found_overlap.get((query_id, target_id), False):
                            find.append(found_overlap[(query_id, target_id)])
                            both.append(found_overlap[(query_id, target_id)])
                        elif found_overlap.get((target_id, query_id), False):
                            find.append(found_overlap[(target_id, query_id)])
                            both.append(found_overlap[(target_id, query_id)])

                    elif all_overlap.get((query_id, target_id), False) or all_overlap.get((target_id, query_id), False):
                        if all_overlap.get((query_id, target_id), False):
                            notfind.append(all_overlap[(query_id, target_id)])
                            both.append(all_overlap[(query_id, target_id)])
                        elif all_overlap.get((target_id, query_id), False):
                            notfind.append(all_overlap[(target_id, query_id)])
                            both.append(all_overlap[(target_id, query_id)])


print("Total pairs: {}, Expected Overlap: {}, Total report: {}, True Positive: {}".format(total_pairs, expect, num_found, true_align))
sensitivity = float(true_align) / expect
accuracy = float(true_align)/num_found
print("sensitivity", sensitivity)
print("FPR", (num_found - true_align) / float(total_pairs- expect))
print("accuracy", accuracy)
print("F1", 2*(accuracy*sensitivity)/(accuracy+sensitivity))
print("Filtration Rate: {}".format(num_found/total_pairs) )

plt.figure()
hist, bins = np.histogram(both, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
hist1, bins = np.histogram(find, bins=bins)
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist1, align='center', width=width, label = "Found", alpha=0.7)
hist1, bins = np.histogram(notfind, bins=bins)
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist1, align='center', width=width, label = "Not found", alpha=0.7)
plt.savefig(fig_path)