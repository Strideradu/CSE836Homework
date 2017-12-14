import numpy as np

c_mat = np.zeros((3, 3), np.int) # first index is predicted, second is real
num_clusters = 0
num_reads = 0

with open("/mnt/home/dunan/Learn/Class/CSE836/virus_reads_k6_32_cluster.out") as f:
    for line in f:
        line = line.strip()
        if line != "" or line[0] != "o":
            count = np.array([0, 0, 0])
            num_clusters += 1

            reads = line.split("\t")
            for read in reads:
                num_reads += 1
                if "HCV" in read:
                    count[0] += 1
                elif "HGV" in read:
                    count[1] += 1
                elif "HXB" in read:
                    count[2] += 1

            max_i = np.argmax(count)
            max_v = count[max_i]

            c_mat[max_i][max_i] += count[max_i]
            for i in range(3):
                if i != max_i:
                    c_mat[max_i][i] += count[i]


print(c_mat)
print("For {} read we generate {} clusters".format(num_reads,num_clusters))
print("{:.2%} HCV read classfied into right clusters (HCV dominated cluster)".format(c_mat[0][0]/np.sum(c_mat[0])))
print("{:.2%} HGV read classfied into right clusters (HGV dominated cluster)".format(c_mat[1][1]/np.sum(c_mat[1])))
print("{:.2%} HIV read classfied into right clusters (HIV dominated cluster)".format(c_mat[2][2]/np.sum(c_mat[2])))
