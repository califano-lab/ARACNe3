# ARACNe3

## ARACNe3 (Algorithm for the Reconstruction of Accurate Cellular Networks ver. 3) in C++

`ARACNe3` is the C++ implementation of `ARACNe`.  `ARACNe3` presents multiple computational and theoretical improvements to the existing `ARACNe-AP` implementation, formulated in Java by Lachmann and colleagues in 2016.  

Lachmann A, Giorgi FM, Lopez G, Califano A. *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information.* **Bioinformatics.** 2016 Jul 15;32(14):2233-5. doi: [10.1093/bioinformatics/btw216](https://dx.doi.org/10.1093/bioinformatics/btw216). Epub 2016 Apr 23.

Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A. *ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context.* **BMC Bioinformatics.** 2006 Mar 20;7 Suppl 1:S7. doi: [10.1186/1471-2105-7-S1-S7](https://dx.doi.org/10.1186/1471-2105-7-S1-S7)

## Download

`git clone https://github.com/califano-lab/ARACNe3-CPP # Clone the repo`

## Building ARACNe3
The C++20 standard is used when compiling `ARACNe3` into an executable.  `ARACNe3` utilizes multiple threads with OpenMP directives, so you must have libraries that support them installed.  Here is _one_ example of how you might download OpenMP libraries on the latest version of MacOS, with homebrew[https://brew.sh] already installed.

`brew install llvm libomp`

Build the executable by cloning the repository and typing `make` in the command line, while in the repository base directory.  You could also manually build the executable by compiling all C++ files and linking `ARACNe3` with object file dependencies in the manner below.  Make sure to utilize all available runtime optimization compiler options.  
### Compiling:
```
g++ -std=c++20 -O3 -c NullModel.cpp; g++ -std=c++20 -O3 -c MatrixReglistIO.cpp; g++ -std=c++20 -O3 -c APMI.cpp; g++ -std=c++20 -O3 -c AlphaPruning.cpp; g++ -std=c++20 -O3 -c MaxEntPruning.cpp; g++ -std=c++20 -O3 -c RegWebFns.cpp; g++ -std=c++20 -O3 -c Consolidator.cpp

```
### Linking:
```
g++ -std=c++20 -O3 ARACNe3.cpp NullModel.o MatrixReglistIO.o APMI.o AlphaPruning.o MaxEntPruning.o RegWebFns.o Consolidator.o -o ARACNe3
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
1.	Normalize a gene expression profile for sequencing depth in each sample.
2.	Run `ARACNe3` according to command line instructions below.

### Optional ways to run ARACNe3
1.	Customizing the population percentage from which to subsample for network generation, or to remove subsampling (`--subsample 1`).
2.	Removing the MaxEnt pruning step will preserve every edge that passes the FDR/FWER Pruning Step.
3.	Customizing the FDR restriction for the first pruning step.  For an FDR restriction, the first pruning step rejects the null hypothesis for edges based on the Benjamini-Hochberg Procedure.  
4. 	Control for the Family Wise Error Rate (FWER) instead of the FDR during the first pruning step (`--FWER`).  If chosen, `--alpha` will be chosen as the FWER control parameter. 

### Output of ARACNe3
`ARACNe3` will output a file for each subnetwork (`-x`) requested with the name "`output_subnet#.txt`" in a directory path provided by the user. Each file shows every significant interaction in three columns:
1.	The regulator.
2.	The target.
3.	The MI (Mutual Information) of the pair.

## Input file format
### A gene/regulator list
A text file, containing one gene symbol per line. Every line has important information. E.g.,
```
g_10_
g_10011_
g_1_
```
### Dataset
A normalized transformed expression profile as a `.tsv` (tab separated value) file, with genes on rows and samples on columns.  Do not include any important information in the first row, except for equal number of columns and spacing as the rows below. E.g.,
```
gene	Sample1	Sample2	Sample3
g_1_	4.99	2.93	0.39
g_10_	0.58	0.18	2.65
g_9432_	3.00	1.27	7.63
g_10006_	1.30	0.05	0.68
g_10011_	0.055	0.73	4.64
```

## Parameters
``-e`` is the expression file

``-r`` is the list of regulators (e.g., TFs)

``-o`` is the output directory

``-x`` is the number of subnetworks to generate (default: `-x 1`)

``-a`` or ``--alpha`` is the alpha parameter for FDR or FWER pruning (default: `--alpha 0.05`)

``--FWER`` tells ARACNe3 to prune by control of FWER instead of the default FDR

``--subsample`` is the population percentage to subsample ($1-e^{-1}$ is default: `--subsample 0.63212...`)

``--noAlpha`` tells ARACNe3 not to prune based on the FDR or FWER (same as: `--alpha 1`)

``--noMaxEnt`` tells ARACNe3 not to prune edges based on the Principle of Maximum Entropy

``--seed`` sets the seed for pseudorandom behavior for null model marginals and subsampling (default: `--seed 0`)

``--noverbose`` removes console messages from ARACNe3 stating elapsed time and resulting edges

## Examples
Note: the examples have been written based on the provided test sets: ``test/exp_mat.txt`` (the normalized expression matrix) and ``test/regulators.txt`` (the list of regulators). 

### Example 1: run ARACNe3 with no pruning steps, no subsampling, and seed equal to 343
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output --subsample 1.00 --noAlpha --noMaxEnt --seed 343
```

### Example 2: run ARACNe3 with all pruning steps, subsampling 33.3% of profiles, controlling for FDR < 0.01
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output --subsample 0.333 --alpha 0.01
``` 

## Currently Under Development:
 - Separating subnetwork generation and consolidation
 - Using the binomial distribution to estimate the significance of an edge, based on how many times it appears in the subnetworks
 - Change entire reg\_reb data structure to use map\_map
 - Use bitwise AND to determine shared edges for MaxEnt pruning step

## Tracking Progress

Whenever a significant change is made to an existing module of this program, test results are noted in the git commit notes.  Standardized comparisons with the Java `ARACNe-AP` are also found in the `test/` directory.  The results reflect runtime performance of the most recent version of `ARACNe3` or a given executable.  Tests were run on a computer with the specifications listed below. 

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

