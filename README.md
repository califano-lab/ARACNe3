# ARACNe3

## ARACNe3 (Algorithm for the Reconstruction of Accurate Cellular Networks ver. 3) in C++

`ARACNe3` is the C++ implementation of `ARACNe`.  `ARACNe3` presents multiple computational and theoretical improvements to the existing `ARACNe-AP` implementation, formulated in Java by Lachmann and colleagues in 2016.  

Lachmann A, Giorgi FM, Lopez G, Califano A. *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information.* **Bioinformatics.** 2016 Jul 15;32(14):2233-5. doi: [10.1093/bioinformatics/btw216](https://dx.doi.org/10.1093/bioinformatics/btw216). Epub 2016 Apr 23.

Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A. *ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context.* **BMC Bioinformatics.** 2006 Mar 20;7 Suppl 1:S7. doi: [10.1186/1471-2105-7-S1-S7](https://dx.doi.org/10.1186/1471-2105-7-S1-S7)

## Download

`git clone https://github.com/califano-lab/ARACNe3-CPP # Clone the repo`

## Building ARACNe3
The C++20 standard is used when compiling `ARACNe3` into an executable.  Build the executable by cloning the repository and typing `make` in the command line (in the repository base directory).  For this, we use the GNU compiler `g++`, so you will need that installed as well.  You could also manually build the executable by compiling all C++ files and linking `ARACNe3` with object file dependencies in the manner below.  Make sure to utilize all available runtime optimization compiler options.  
### Compiling:
```
g++ -std=c++20 -O3 -c NullModel.cpp; g++ -std=c++20 -O3 -c MatrixReglistIO.cpp; g++ -std=c++20 -O3 -c APMI.cpp; g++ -std=c++20 -O3 -c FDRPruning.cpp; g++ -std=c++20 -O3 -c MaxEntPruning.cpp; g++ -std=c++20 -O3 -c RegWebFns.cpp 

```
### Linking:
```
g++ -std=c++20 -O3 ARACNe3.cpp NullModel.o MatrixReglistIO.o APMI.o FDRPruning.o MaxEntPruning.o RegWebFns.o -o ARACNe3
``` 
### Troubleshooting Build Issues
If you face issues building `ARACNe3` besides lacking a C++20 compiler, or if you are experiencing `segmentation fault` even with proper usage, try compiling `ARACNe3` without the GNU `-O3` compiler optimization option, or repace it with `-O2`. 

## Using ARACNe3
### Input files needed to run ARACNe3
See below for file format specification (or download the test files from our repository)
1.	A normalized expression matrix (CPM, TPM, etc.)
2.	List of regulators (e.g. Transcription Factors)

*Note: All regulators must have a defined expression profile in the expression matrix*.

### Steps required to run ARACNe3
1.	Normalize a gene expression profile for sequencing depth in each sample, which could be CPM or TPM
2. 	Run `ARACNe3` according to command line instructions below

### Optional ways to run ARACNe3
1.	Customizing the population percentage from which to subsample for network generation, or to remove subsampling (`--subsample 1.00`).
1.	Removing the MaxEnt pruning step will preserve every edge that passes the FDR Pruning Step
2.	Customizing the FDR restriction for the first pruning step, which rejects the null hypothesis for edges based on the Benjamini-Hochberg Procedure.

### Output of ARACNe3
`ARACNe3` will output a file `output_(FINAL NETWORK SIZE).txt` in a directory path provided by the user (e.g. `test/output/output_127304.txt`). This file shows every significant interaction in three columns:
1.	The regulator.
2.	The target.
3.	The MI (Mutual Information) of the pair.

## Input file format
### A gene/regulator list
A text file, containing one gene symbol per line. Every line has important information. E.g.,
```
g_9970_
g_9971_
g_9975_
g_9984_
g_9987_
```
### Dataset
A normalized transformed expression profile as a `.tsv` (tab separated value) file, with genes on rows and samples on columns.  Do not include any important information in the first row, except for equal number of columns and spacing as the rows below. E.g.,
```
gene    Sample1   Sample2   Sample3
g_1_	4.99	2.93	0.39
g_10_   0.58       0.18       2.65       0.73
g_10006_        1.30        0.05      0.68
g_10011_        0.055      0.73       4.64
```

## Parameters
``-e`` is the expression file

``-r`` is the list of regulators (e.g., TFs)

``-o`` is the output directory

``--subsample`` is the population percentage to subsample (default `--subsample 0.6321205588`)

``--FDR`` is the FDR parameter to set (default: `--FDR 0.05`)

``--noFDR`` tells ARACNe3 not to prune based on the FDR (same as: `--FDR 1.00`)

``--noMaxEnt`` tells ARACNe3 not to prune edges based on the Principle of Maximum Entropy

``--noverbose`` removes console messages from ARACNe3 stating elapsed time and resulting edges

## Examples
Note: the examples have been written based on the provided test sets: ``test/exp_mat.txt`` (the normalized expression matrix) and ``test/regulators.txt`` (the list of regulators). 

### Example 1: run ARACNe3 with no pruning steps, no subsampling
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output --subsample 1.00 --noFDR --noMaxEnt
```

### Example 2: run ARACNe3 with all pruning steps, subsampling 33.3% of profiles, controlling for FDR < 0.01
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output --subsample 0.333 --FDR 0.01
``` 

## Currently Under Development:
 - Replace `test/exp_mat.txt` with a normalized expression profile, as opposed to a copula-transformed one
 - Multithreading/non-multithreading option using standard library
 - Remove any references on data types `<=4B`, as references instantiate pointer values which are at least 4B (and typically are 8B on 64-bit systems)
 - Return by reference when applicable!
 - Optimize p-value calculation for each MI value
 - Low-level optimization and parallel for loop processing. Namely, minimizing heap allocation and using caches, as well as using the most efficient data structures required to store edge information (hashmaps, linked lists, adjacency matrices, etc.)
 - Can using prime products of each regulator target set immediately identify which targets are shared?  Explore making O(N^3) MaxEnt pruning step into an O(N^2) step
 
 Plan of action: Overhaul of Low-level algorithms -> Rewriting codebase with highest efficiency class design and matured codebase -> Implement multithreading 

## Tracking Progress

Whenever a significant change is made to an existing module of this program, test results are noted in the git commit notes.  Standardized comparisons with the Java `ARACNe-AP` are also found in the `test/` directory.  The results reflect runtime performance of the most recent version of `ARACNe3` or a given executable.  Tests were run on a computer with the specifications listed below. 

## System Information
| System Feature | Value |
| :----: | :----: |
| Model Identifier: | MacBookPro18,1 |
| Chip: | Apple M1 Pro |
| Total Number of Cores: | 10 (8 performance and 2 efficiency) |
| Memory: | 32 GB |
| System Firmware Version: |7459.101.3 |
| OS Loader Version: | 7459.101.3 |

## Contact
Please contact Aaron Griffin, Lukas Vlahos, or Andrew Howe for questions regarding this project.

Aaron T. Griffin - atg2142@cumc.columbia.edu 

Lukas J. Vlahos - lv2395@cumc.columbia.edu 

Andrew R. Howe - arh2207@columbia.edu

Previous commits and developments can be found in the following repository:
```
https://github.com/arhowe00/ARACNe3_Prototypes
```

