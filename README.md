# PFP-eBWT

PFP-eBWT is a tool that computes the eBWT of string collections using the Prefix-Free Parsing (PFP) method.

Building upon the original implementation, this repository introduces additional functionality:

- Computation of the document array (DA) associated with the eBWT.
- A compact data structure that enables eBWT inversion while preserving the original order of strings in the collection.
- A script to compute different dissimilarity measures pairwise on a string collection, leveraging the document array and the eBWT, and to build the corresponding distance matrix.

## Installation

### Download and Build

```bash
git clone https://github.com/SaraCaporale/PFP-eBWT.git
cd PFP-eBWT
mkdir build
cd build
cmake ..
make
g++ ../distances.cpp -o distances
```

### Run on Example Data

```bash
# Build the eBWT on a toy data set
python3 pfpebwt ../yeast.fasta -w 10 -p 100

# Calculate dissimilarity measure on toy dataset
./distances -r ../yeast.fasta ../yeast.fasta.da ../yeast.fasta.ebwt
```

## Usage

### Construction of the eBWT and the corresponding DA:

```
usage: pfpebwt [-h] [-w WSIZE] [-p MOD] [-t T] [-n N] [-r] [-v] [-d] [--reads]
               [--period] [--invert] [--keep] [--parsing]
               input

positional arguments:
  input                 input fasta file name

optional arguments:
  -h, --help            show this help message and exit
  -w WSIZE, --wsize WSIZE
                        sliding window size (def. 10)
  -p MOD, --mod MOD     hash modulus (def. 100)
  -t T                  number of helper threads (def. None)
  -n N                  number of different primes (def. 2)
  -r                    store the eBWT in RLE format (def. False)
  -v                    verbose (def. False)
  -d                    use different remainders instead of different prime numbers (def. False)
  --reads               process input as a read collection (def. False)
  --period              remove sequences which are not primitive (def. False)
  --invert              invert the eBWT (def. False)
  --keep                keep auxiliary files (debug only)
  --parsing             stop after the parsing phase (debug only)
```

The default algorithm will run with one prime number and one remainder to search for the trigger strings. If you want to build the eBWT of a collection of short sequences you have to use the `--reads` flag, and set with `-n` the maximum number of different primes allowed. With `-d` you can use different remainders instead of different primes. 
You can activate the `--period` flag in case your dataset contains non-primitive words, this flag will allow to filter them out. 
You can activate the `--invert` flag in case you want to invert the eBWT and write the inverted sequences to a file, this process requires to compute the WT of the eBWT in the internal memory. 

The computation of the document array associated with the eBWT of the string collection has been integrated into the default algorithm. In addition to the files created by the original implementation, this version generates the following additional files with extension:
- `.da`: contains the DA associated with the eBWT of the string collection,
- `.dap`: contains the DA related to the eBWT of the parsing, used as an auxiliary data structure for computating the DA.

The original inversion algorithm returns the string collection in alphabetic order This version implents a method to invert the eBWT while preserving the original order of the strings in the collection. The algorithm generates a file with extension:
- `.bstOrder`: contains the permutation of the indexes of the original strings in the eBWT. 

### Compute dissimilarity measures

```
usage: ./distances <-d | -r | -m | -e> labels_filename DA_filename [ebwt_filename (only with -r)]

help:
The following flags are **mutually exclusive**
  -d: computes the delta measure [6]
  -r: computes the rho measure [5]
  -m: computes the d_M measure [7] 
  -e: computes the d_E measure [7]

parameters:
  labels_filename:    
    file that contains the labels that identify the strings in
    the collection. The file must follow a similar format to the 
    fasta format, i.e. each label should be on a different line 
    and preceded by the symbol '>', all lines that don't start 
    with the symbol '>' will be ignored. In fact, labels_filename 
    could be the name of the fasta file given as input to the 
    pfpebwt algorithm. However, if the fasta file is big, 
    consider creating a file containing only the labels (in the 
    correct order), so as to speed up the computation.
  DA_filename:
    file containing the DA associated with the eBWT of the string 
    collection
  ebwt_filename:
    (to use only if -r is set) file containing the eBWT of the 
    string collection 
```

The algorithm will compute the dissimilarity measure specified by the flag and will create a file containing the corresponding distance matrix in PHYLIP format.


## External resources

* [gSACA-K](https://github.com/felipelouza/gsa-is.git)
* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)

## Authors

### Theoretical results:

* Christina Boucher
* Davide Cenzato
* Zsuzsanna Lipták
* Massimiliano Rossi
* Marinella Sciortino

### Implementation and experiments:

* [Davide Cenzato](https://github.com/davidecenzato) 
* [Massimiliano Rossi](https://github.com/maxrossi91)

# References

[1] Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini: *Prefix-Free Parsing for Building Big BWTs.* In Proc. of the 18th International Workshop on Algorithms in Bioinformatics, WABI 2018.

[2] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini and Taher Mun: *Prefix-free parsing for building big BWTs.* Algorithms Mol. Biol. 14(1): 13:1-13:15 (2019)

[3] Ge Nong, Sen Zhang and Wai Hong Chan: *Two Efficient Algorithms for Linear Time Suffix Array Construction.* IEEE Trans. Computers 60(10): 1471-1484 (2011)

[4] Hideo Bannai, Juha Kärkkäinen, Dominik Köppl and Marcin Piatkowski: *Constructing the bijective and the extended Burrows-Wheeler-Transform in linear time.* In Proc. of the 32nd Annual Symposium on Combinatorial Pattern Matching, CPM 2021.

[5] S Mantaci et al. “A New Combinatorial Approach to Sequence Comparison”. In: Theory of Computing Systems 42.3 (apr. 2008), pp. 411–429.

[6] Sabrina Mantaci et al. “An Extension of the Burrows Wheeler Transform and Applications to Sequence Comparison and Data Compression”. In: Combinatorial Pattern Matching. A cura di Alberto Apostolico, Maxime Crochemore e Kunsoo Park. Berlin, Heidelberg: Springer Berlin Heidelberg, 2005, pp. 178–189. isbn: 978-3-540-31562-9.

[7] Lianping Yang, Xiangde Zhang e Tianming Wang. “The Burrows–Wheeler similarity distribution between biological sequences based on Burrows–Wheeler transform”. In: Journal of Theoretical Biology 262.4 (2010), pp. 742–749. issn: 0022-5193. doi: https://doi.org/10.1016/j.jtbi.2009.10.033. url: https://www.sciencedirect.com/science/article/pii/S0022519309005220.

[8] Christina Boucher et al. “Computing the Original eBWT Faster, Simpler, and with Less Memory”. In: Proceedings of 28th International Symposium in String Processing and Information Retrieval SPIRE 2021. Vol. 12944. Lecture Notes in Computer Science. 2021, pp. 129–142.
