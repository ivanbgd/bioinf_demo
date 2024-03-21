# A Bioinformatics Demo in Python

## Description

### The Goal
The goal was to write a new program without using existing tools to solve the following problem.  
We have a FASTQ file in the GZIP format and an adapter sequence list as the inputs.  
We would like to implement some filtering on the FASTQ file and output clean results also in the GZIP format.  
The filtering steps include poly (X) detection and adapter search.  
Aside of FASTQ, we'd also like to output the number of each filter item, such as how many reads were filtered by polyX 
or by an adapter.  
We wish to have a high-performance design based on multiple threads for I/O and for the algorithm of sequences search.

### The Filtering Process (in order)
1. Search for polyX whose length is less than or equal to 15 in the sequence, allowing 1 mismatch.
If poly-X is found, the read will be discarded.
2. If the poly-X is not found, then search for the adapter whose length is less than or equal to 32 bp
in the sequence without any mismatch.  
If the adapter is found, which means it is identical to a part of the sequence, the record is to be discarded.
3. If the adapter is not found either, the read shall be considered as a clean one and shall be output to the
result file.
4. Count how many records are filtered due to the existence of poly-X and adapter respectively.

### Test data
- Inputs: a FASTQ file (input.fq.gz) and an adapter sequence list (adapter.list)  
- Outputs: a FASTQ file (out.fq.gz) and filter stat (out.stat.txt)  

### Glossary
- [FASTQ](http://maq.sourceforge.net/fastq.shtml): Stores sequences and Phred qualities in a single file.  
- polyX: Include polyA, polyC, polyG, polyT, for example: AAAAAAAAAAAAAAA, CCCCCCCCCCCCCC.
- Adapter: A short DNA sequence that includes A/C/G/T bases which should be trimmed or discarded from a read.

## Activating Virtual Environment
On Windows:  
`.\venv\Scripts\activate`

## Running the program
Run the Python interpreter and the program from the main project directory, named *bioinf_demo*.  
The main "executable" file is located inside the *bin* subdirectory.  
`python bin/bioinf_demo.py`

## Implementation
- Has both Naive and Prefix Trie Matching
- Has both sequential and parallel implementations
  - [*Modin*](https://modin.readthedocs.io/en/stable/)

Modin is **very** easy to translate to from Pandas.
It preserves the order of records in the output file. 
Modin can be used on a single computer (shared-memory model), or in a cluster (distributed-memory model).
It can read and process larger-than-memory files.

## Remarks & Notes
- FASTQ files are mostly used to store short-read data from high-throughput sequencing experiments. 
  The sequence and quality scores are usually put into a single line each, and indeed many tools assume that
  each record in a FASTQ file is exactly four lines long, even though this isn’t guaranteed.
  [Source](https://bioinformatics.stackexchange.com/questions/14/what-is-the-difference-between-fasta-fastq-and-sam-file-formats).

## Potential Improvements
- [*biopython*](https://biopython.org/) supports [*PyPy*](https://www.pypy.org/) as documented on its
  [PyPI](https://pypi.org/project/biopython/) page.
- [*pyfastx*](https://pyfastx.readthedocs.io/en/latest/) didn't support Python 3.10 at the time of experimenting;
  it didn't want to install.

## Run Tests

### Unit Tests

The `unittest` runner will run both unit tests and doctests if they exist.

`python -m unittest`

For verbose output:

`python -m unittest -v`
