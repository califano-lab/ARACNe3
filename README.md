# ARACNe3

## ARACNe3 (Algorithm for the Reconstruction of Accurate Cellular Networks ver. 3)
ARACNe3 is an implementation of ARACNe which presents computational improvements and theoretical changes to the recent ARACNe-AP implementation.  Given a list of regulators and a gene expression profile, ARACNe3 is used to infer irreducibly dependent regulatory interactions and output a Gene Regulatory Network (GRN). The mainstream analysis generates many GRNs of subsamples of the profile ("subnetworks") and then consolidates the subnetworks into a robust GRN.  The consolidated output is a GRN whose edge strengths can be quantified by several statistical metrics.

Lachmann A, Giorgi FM, Lopez G, Califano A. *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information.* **Bioinformatics.** 2016 Jul 15;32(14):2233-5. doi: [10.1093/bioinformatics/btw216](https://dx.doi.org/10.1093/bioinformatics/btw216). Epub 2016 Apr 23.

Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A. *ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context.* **BMC Bioinformatics.** 2006 Mar 20;7 Suppl 1:S7. doi: [10.1186/1471-2105-7-S1-S7](https://dx.doi.org/10.1186/1471-2105-7-S1-S7)

## Downloading ARACNe3
You may [download ARACNe3](https://github.com/califano-lab/ARACNe3) for MacOS and Windows or build ARACNe3 manually for niche uses, like on an HPC cluster.

If running the ARACNe3 executables fails, on MacOS try installing the [Xcode Command Line Tools](https://mac.install.guide/commandlinetools/4.html), and on Windows install the [Visual C++ Redistributable for Visual Studio](https://docs.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist).  

## Building ARACNe3 (Optional)
### Installing Libraries
ARACNe3 supports multithreading using OpenMP.  Here is **one** example of how you might download OpenMP libraries on the latest version of MacOS, with [homebrew](https://brew.sh) already installed.

```
brew install libomp  # Install OpenMP
```

### Build ARACNe3 with CMake
`cmake` is used to build ARACNe3 and simplify cross-platform compatibility. After installing `cmake`, follow the instructions below to build ARACNe3.

```
git clone https://github.com/califano-lab/ARACNe3  # Clone the repo

cd ARACNe3/
git submodule update --init --recursive  # Install dependencies

mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build  # Build ARACNe3
```

The executable should be built in `./build/src/app/ARACNe3_app_release`. If you are having issues, try using `cmake3` instead of `cmake` in the instructions above.

## Using ARACNe3
### Steps required to run ARACNe3
1.	Normalize a gene expression profile for sequencing depth in each sample (CPM, TPM, etc.).
2.	Run ARACNe3 according to command line instructions below.

### ARACNe3 input files (see [input file format](#input-file-format) below)
1.	A `tsv` that contains the `G+1 x N+1` normalized expression profile, with genes as rows and samples as columns. Note: the "`+1`" is for header names & column names, respectively.
2.	List of regulators (e.g., transcription factors).

### ARACNe3 output files
ARACNe3 outputs in the directory provided by the user (e.g., `-o outputdir/`).  If ARACNe3 is run with `--runid abc`, within `outputdir/` the subnetworks are in `outputdir/subnets/`, the consolidated network is the file `outputdir/consolidated-net_abc.tsv`, and log information for each subnetwork is in `outputdir/subnets_log/`. Log information for the overall ARACNe3 instance is the file `outputdir/log_abc.txt`.  

The subnetworks directory `outputdir/subnets/` contains a file for each subnetwork named `subnet#_abc.tsv`.  Each subnetwork file describes significant interactions in three columns:
1.	The regulator.
2.	The target.
3.	The APMI (Mutual Information estimated by Adaptive Partitioning) of the pair, based on the subsampled gene expression profiles.

The consensus network `outputdir/consolidated-net_abc.tsv` is a file that consolidates subnetworks and summarizes each edge in five columns:
1.	The regulator.
2.	The target.
3.	The APMI of the pair, based on the full gene expression profiles.
4.	The Spearman's Rank Correlation Coefficient of the pair, based on the full gene expression profile.
5.	The _p_-value of the edge, based on how many subnetworks in which it was identified (see the [manuscript](https://github.com/califano-lab/ARACNe3) for methodology).

## Parameters
### Required
`-e` is the expression file.

`-r` is the list of regulators (e.g., TFs).

`-o` is the output directory.

### Useful
`-x` adjusts the stopping criteria when generating or using multiple subnetworks.  By default (`-x 1`) specifies the fixed number of subnetworks to generate, in this case just 1.

`--adaptive` changes the meaning of `-x` to specify regulon occupancy, instead of number of subnetworks to generate.  Regulon occupancy is defined as the minimum number of unique targets observed per regulator, when all subnetworks are consolidated into one (default: `-x 30` if `--adaptive` is specified).  Note that different subnetworks find different targets because of different subsamples, but in the consolidated network there are metrics for confidence in these targets.

`--alpha` is the cutoff for False Discovery Rate (FDR) pruning (default: `--alpha 0.05`).  The FDR cutoff is the first pruning step, which rejects the null hypothesis for edges based on the Benjamini-Hochberg Procedure.

`--seed` fixes the seed for deterministic behavior (e.g.: `--seed 9001`).

`--threads` sets the number of threads to use (default: `--threads 1`).

`--runid` allows you to pass an identifier to replace `defaultid` in `log_defaultid.txt`. Does not affect each modulator's log, only the instance log (default: `--runid defaultid`).

### Optional
`--FWER` tells ARACNe3 to prune by control of Family Wise Error Rate (FWER) alpha, instead of the default control for FDR.

`--FPR` tells ARACNe3 to prune by control of False Positive Rate (FPR) alpha (a.k.a., _p_-value), instead of the default control for FDR.

`--subsample` is the population percentage to subsample ($1-e^{-1}$ is default: `--subsample 0.63212`).

`--noAlpha` tells ARACNe3 to skip subnetwork pruning by FDR/FPR/FWER (same as: `--alpha 1`).

`--noMaxEnt` tells ARACNe3 to skip subnetwork pruning by Principle of Maximum Entropy, the second pruning step.

`--noConsolidate` tells ARACNe3 not to consolidate subnetworks, only keeping the final log, the `log/` subdirectory, and all subnetworks generated in `subnets/`.

`--consolidate` tells ARACNe3 to skip generating subnetworks and consolidate existing subnetworks.  An expression file and a list of regulators must still be provided with `-e` and `-r`, respectively.  `-o` specifies the directory location of an ARACNe3 output.  Finally, `-x` specifies how many subnetwork files to use in consolidate (default: `-x 1`). Note that output directory `-o` _**must**_ contain the subdirectories `subnets/` and `log/` that follow the exact conventions as an ARACNe3 output (including numbering).  Each subnetwork used must be mapped 1:1 with its log file because consolidation generates _p_-values for edges strictly based on parameters used during the subnetwork generation, which are stored in the log files.

## Examples
Note: the examples have been written based on the provided test sets: `test/exp_mat.txt` (the normalized expression matrix) and `test/regulators.txt` (the list of regulators).

### Example 1: generate one ARACNe3 subnetwork, subsampling $1-e^{-1}%$ of expression profiles, controlling for FDR < 0.05, using only one CPU core
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output
```

### Example 2: generate one ARACNe3 subnetwork with no pruning steps, no subsampling, and seed equal to 343, using 10 CPU cores
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output --subsample 1.00 --noAlpha --noMaxEnt --seed 343 --threads 10
```

### Example 3: generate one ARACNe3 subnetwork with all pruning steps, subsampling 33.3% of profiles, controlling for FDR < 0.01
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output --subsample 0.333 --alpha 0.01
``` 

### Example 4: generate thirty ARACNe3 subnetworks, subsampling $1-e^{-1}%$ of expression profiles, controlling for FDR < 0.05
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output -x 30
``` 

### Example 5: generate ARACNe3 subnetworks adaptively, until at least 50 targets are observed per regulon in the consensus network, controlling for FWER < 0.10, and skipping the MaxEnt pruning step
```
./ARACNe3 -e test/exp_mat.txt -r test/regulators.txt -o test/output -x 50 -FWER --alpha 0.10 --adaptive --noMaxEnt
``` 

## Input file format
### A gene/regulator list
A text file, containing one symbol per line. E.g.,

```
g_10_
g_10011_
g_1_
```

### Expression file
A `G+1 x N+1`, normalized expression profile, with genes on rows and samples on columns (Note: the `+1` is extra, from the row names and column names). This should be in `tsv` format and have both row and column names. Row names are essential for ARACNe3 to store a gene's expression. Column names are used to count the number of samples. For example, the `5+1 x 3+1` matrix below is a compliant `tsv`:

```
col	names	NOT	important
g_1_	4.99	2.93	0.39
g_10_	0.58	0.18	2.65
g_9432_	3.00	1.27	7.63
g_10006_	1.30	0.05	0.68
g_10011_	0.055	0.73	4.64
```

## Contact
Please contact Aaron Griffin (theory) or Andrew Howe (codebase) for questions regarding this project.

Aaron T. Griffin - atg2142@cumc.columbia.edu

Andrew R. Howe - arh2207@columbia.edu
