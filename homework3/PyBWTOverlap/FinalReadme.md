###PyBWTOverlap

#### To Build Index

* I use python3 for this homework. To run the program, use following command

    
    usage: PyBWTBuild.py [-h] input save

    positional arguments:
      input       path for the fasta file to build bwt
      save        path to save bwt
    
    optional arguments:
      -h, --help  show this help message and exit

#### To Recruit All Overlap Reads
* Run follow command, and the program will generate a fasta file contain all reads it found has overlap:

    
    usage: PyBWTFinaAll.py [-h] seed data index threshold output
    
    positional arguments:
      seed        path for query sequence pasta file
      data        path for all reads
      index       path for bwt index
      threshold   threshold to report overlap
      output      path for output fasta
    
    optional arguments:
      -h, --help  show this help message and exit



#### Output Format For Overlap
The output of program has following fields and using tab delimited

    query id, query start, query end, target id, target start, target end, overlap size, direction of query