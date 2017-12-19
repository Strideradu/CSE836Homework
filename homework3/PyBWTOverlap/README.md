###PyBWTOverlap

#### To Run

* I use python3 for this homework. To run the program, use following command

    
    python PyBWOverlap.py query target threshold
    
    positional arguments:
        query       path for query sequence pasta file
        target      path for target sequence pasta file
        threshold   threshold to report overlap
    
    optional arguments:
        -h, --help  show this help message and exit

* To use python3 on arctics.cse.msu.edu:

    
    python3 PyBWOverlap.py query target threshold

* On MSU HPCC, load python 3 module first, And then use the python command
    
    
    module load Python/3.5.3
    python PyBWOverlap.py query target threshold



query id, query start, query end, target id, target start, target end, overlap size, direction of query#### Output Format
The output of program has following fields and using tab delimited

    
