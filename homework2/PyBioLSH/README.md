###PyBioLSH

#### To Run

I use python3 for this homework. To run the program, use following command (in arctics you should replace python with python3)

    python PyBioLSH.py [-h] [-cluster] path k perm threshold
    
    positional arguments:
      path        path for sequence pasta file
      k           size of kmer
      perm        number of permutation
      threshold   threshold to report similar
    
    optional arguments:
      -h, --help  show this help message and exit
      -cluster    generate clusters instead of find similar reads
      
if ```-cluster ``` used, then the threshold parameters will be ignored

The default hashing function is refered from https://naml.us/post/inverse-of-a-hash-function/
The LSH and MinHash modified from https://github.com/ekzhu/datasketch
And the fasta/fastq reader is from https://github.com/lh3/readfq