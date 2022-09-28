# ARACNe3

## ARACNe3 (Algorithm for the Reconstruction of Accurate Cellular Networks ver. 3) in C++

`ARACNe3` is the C++ implementation of `ARACNe`.  `ARACNe3` presents multiple computational and theoretical improvements to the existing `ARACNe-AP` implementation, formulated in Java by Lachmann and colleagues in 2016.  

Lachmann A, Giorgi FM, Lopez G, Califano A. *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information.* **Bioinformatics.** 2016 Jul 15;32(14):2233-5. doi: [10.1093/bioinformatics/btw216](https://dx.doi.org/10.1093/bioinformatics/btw216). Epub 2016 Apr 23.

Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A. *ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context.* **BMC Bioinformatics.** 2006 Mar 20;7 Suppl 1:S7. doi: [10.1186/1471-2105-7-S1-S7](https://dx.doi.org/10.1186/1471-2105-7-S1-S7)

## Downloading ARACNe3

`git clone https://github.com/califano-lab/ARACNe3-CPP # Clone the repo`

After cloning the repository, you may build `ARACNe3` manually or access the latest executables in the `products/` directory.  Currently, `ARACNe3` has been pre-compiled to work on MacOS 12.0 or later (in `products/macOS 12.0/`) and Windows 10 (`products/windows/`).  Make sure that if using Windows, you choose the version compiled for your CPU architecture (ARM 64-bit or x86 64-bit)

Note that you may require additional software to run the C++ executables.  If running the `ARACNe3` executables fails, on MacOS try installing the [Xcode Command Line Tools](https://mac.install.guide/commandlinetools/4.html), and on Windows install the [Visual C++ Redistributable for Visual Studio](https://docs.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist).  

## Building ARACNe3 Manually (Optional)
The C++17 standard is used when compiling `ARACNe3` into an executable.  `ARACNe3` utilizes multiple threads with `OpenMP` directives, so you must have libraries that support them installed.  Here is **one** example of how you might download OpenMP libraries on the latest version of MacOS, with [homebrew](https://brew.sh) installed already.

`brew install llvm libomp`

Build the executable by cloning the repository, compiling all C++ files, and then linking `ARACNe3` with object file dependencies in the manner below.  Make sure that `OpenMP` libraries have been installed, and their include directories and library directories have been implied by your [environment variables](https://en.wikipedia.org/wiki/Environment_variable).  Make sure to utilize all available runtime optimization compiler options.  
### Compiling Example:
```
g++ -std=c++17 -O3 -fopenmp -c ARACNe3.cpp NullModel.cpp IO.cpp APMI.cpp AlphaPruning.cpp MaxEntPruning.cpp RegWebFns.cpp Consolidator.cpp -lstdc++fs
```
### Linking Example:
```
g++ -std=c++17 -O3 -fopenmp ARACNe3.o NullModel.o IO.o APMI.o AlphaPruning.o MaxEntPruning.o RegWebFns.o Consolidator.o -o ARACNe3 -lstdc++fs
``` 
### Troubleshooting Build Issues
If you face issues building `ARACNe3` besides lacking a C++17 compiler, or if you are experiencing `segmentation fault` even with proper usage, try compiling `ARACNe3` without the clang `-O3` compiler optimization option, or repace it with `-O2`. 

## Using ARACNe3
### Input files needed to run ARACNe3
See below for file format specification (or download the test files from our repository)
1.	A `G+1 x N+1` normalized expression matrix (CPM, TPM, etc.)
2.	List of regulators (e.g. Transcription Factors)

### Steps required to run ARACNe3
1.	Normalize a gene expression profile for sequencing depth in each sample.
2.	Run `ARACNe3` according to command line instructions below.

### Optional ways to run ARACNe3
1.	Customizing the population percentage from which to subsample for network generation, or to remove subsampling (`--subsample 1`).
2.	Removing the MaxEnt pruning step will preserve every edge that passes the FDR/FWER Pruning Step.
3.	Customizing the FDR restriction for the first pruning step.  For an FDR restriction, the first pruning step rejects the null hypothesis for edges based on the Benjamini-Hochberg Procedure.  
4.	Control for the Family Wise Error Rate (FWER) instead of the FDR during the first pruning step (`--FWER`).  If chosen, `--alpha` will be chosen as the FWER control parameter. 
5.	Using regulon occupancy, or minimum targets per regulator, as a stopping criteria as opposed to number of subnetworks (`--adaptive -x 30`).

### Output of ARACNe3
`ARACNe3` will output a directory structure within the output directory provided by the user (e.g., `-o outputdir/`) that contains the subnetworks (in `outputdir/subnets/`), the consensus network (`outputdir/finalNet.txt`), and log information for the subnetworks (in `outputdir/log/`) and the overall `ARACNe3` runtime (`outputdir/finalLog.txt`).  The subnetworks subdirectory will contain a file for each subnetwork (`-x`) requested, automatically named `output_subnet#.txt`. Each file shows every significant interaction in three columns:
1.	The regulator.
2.	The target.
3.	The APMI (Mutual Information estimated by Adaptive Partitioning) of the pair, based on the subsampled gene expression profiles of the regulator and target.

The consensus network `outputdir/finalNet.txt` is a file that consolidates all generated networks and summarizes each edge in five columns:
1.	The regulator.
2.	The target.
3.	The APMI of the pair, based on the full gene expression profiles of the regulator and target.
4.	The Spearman's Rank Correlation Coefficient of the pair, based on the full gene expression profile.
5.	The p-value of the edge, based on how many subnetworks in which it was identified (see the [manuscript](https://github.com/califano-lab/ARACNe3-CPP) for methodology).

## Input file format
### A gene/regulator list
A text file, containing one gene symbol per line. Every line has important information. E.g.,
```
g_10_
g_10011_
g_1_
```
### Dataset
A `G+1 x N+1`, normalized and transformed expression profile, with genes on rows and samples on columns.  All columns should be tab-delimited (`\t`), and all rows should be newline delimited (`\n`).  Do not include any important information in the first row, except for an equal number of columns (delimited by tab) as the rows below. E.g.,
```
gene	Samp5	Sample2	Samp1
g_1_	4.99	2.93	0.39
g_10_	0.58	0.18	2.65
g_9432_	3.00	1.27	7.63
g_10006_	1.30	0.05	0.68
g_10011_	0.055	0.73	4.64
```

## Parameters
``-e`` is the expression file.

``-r`` is the list of regulators (e.g., TFs).

``-o`` is the output directory.

``-x`` is the stopping criteria for generating multiple subnetworks or consolidating subnetworks.  By default, it specifies a fixed number of subnetworks to generate (default: `-x 1`).

``--adaptive`` changes the stopping criteria `-x` to specify regulon occupancy, instead of number of subnetworks to generate.  Regulon occupancy is defined as the minimum number of unique targets observed per regulator, if all subnetworks are consolidated into one (default: `-x 30`).

``--alpha`` is the alpha parameter for FDR or FWER pruning (default: `--alpha 0.05`).

``--FWER`` tells ARACNe3 to prune by control of FWER alpha, instead of the default control for FDR alpha.

``--subsample`` is the population percentage to subsample ($1-e^{-1}$ is default: `--subsample 0.63212...`).

``--noAlpha`` tells ARACNe3 not to prune based on the FDR or FWER (same as: `--alpha 1`).

``--noMaxEnt`` tells ARACNe3 not to prune edges based on the Principle of Maximum Entropy.

``--seed`` sets the seed for all programmatic pseudorandom behavior (default: `--seed 0`).

``--threads`` sets the number of threads to use during subnetwork generation (default: `--threads 1`).

``--noconsolidate`` tells ARACNe3 to skip the consolidation step and instead only preserve the null model, the final log, the `log/` subdirectory, and all subnetworks generated in `subnets/`.

``--consolidate`` tells ARACNe3 to run in consolidate mode.  An expression file and a list of regulators still must be provided with `-e` and `-r`, respectively.  `-o` specifies the directory location.  Finally, `-x` specifies how many subnetwork files to use (default: `-x 1`). Note that output directory structure `-o` **must** contain subdirectories `subnets/` and `log/` that follow the exact conventions as an ARACNe3 output.

## Examples
Note: the examples have been written based on the provided test sets: ``test/exp_mat.txt`` (the normalized expression matrix) and ``test/regulators.txt`` (the list of regulators).  Due to the limited capacity of the repository, you may need to decompress the expression matrix before usage.

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

## Contact
Please contact Aaron Griffin for questions regarding this project.

Aaron T. Griffin - atg2142@cumc.columbia.edu 

Previous commits and developments can be found in the following repository:
```
https://github.com/arhowe00/ARACNe3_Prototypes
```
