# ARACNe3 C++

## ARACNe3 (Algorithm for the Reconstruction of Accurate Cellular Networks ver. 3) in C++

This is a prototype of the C++ implementation of ARACNe3 currently under development.  ARACNe3 presents multiple computational and theoretical improvements to the existing ARACNe-AP algorithm formulated in Java.  


Please contact Aaron T. Griffin and Lukas J. Vlahos for any questions regarding this project.


Aaron T. Griffin - atg2142@cumc.columbia.edu 

Lukas J. Vlahos - lv2395@cumc.columbia.edu 

## Running the Most Recent Executable

This codebase is currently under development, but the most up-to-date ARACNe3 executable can still be run on any platform.  The executable must be compiled using source files in this directory, but only standard libraries are used, so this can be done by linking all object file dependencies with the ARACNe3 target executable, which can be run using the following commandline arguments: 


`./ARACNe3 /path/to/regulators.txt /path/to/gexpmatrix.txt`

The program will output an 'output.txt' that contains the regulator-target MI values estimated via Adaptive Partitioning (APMI), and the associated p-value determined from a null distribution of the APMI between 1 million shuffled gene-expression marginals.  Currently, this is computationally slow, can only run on one CPU core, and will take approximately 12 minutes on a computer with the specifications below.   

## List of Improvements in Development:
 - **IN PROGRESS** First MI-pruning step based on null model for MI and Benjamini-Hochberg principle control for user-parametrizable FDR (p values done, must make function to vectorize the p-values and sort, then prune)
 - Second DPI-pruning step made optional and implemented
 - Optimize calculation of 1 million null MIs, possibly by shuffling index and sending to APMI(), as opposed to vector matrix by std::shuffle()
 - Multithreading/non-multithreading option using standard library
 - Compression of gene name identifiers ("_gxxx_") from std::string -> uint16_t and decompression
 - Optimize p-value calculation for each MI value
 - Low-level optimization and parallel for loop processing. Namely, minimizing heap allocation and using caches, as well as using the most efficient data structures required to store edge information (hashmaps, linked lists, adjacency matrices, etc.)

## Tracking Progress

Whenever a significant change is made to an existing module of this program,
such as MatrixReglistIO.cpp for reading tsv or regulator lists and forming data
structures, or NullModel.cpp for creating the null distribution for mutual information, test results are appended to test/current\_test.txt.  The results
reflect the command run on the most recent version of the given executable, such
as ./NullModel for testing computation of the null model for MI.  Tests were run
on the following computer, and the 'time' program was released in MacOS 12.3 as
a utility expected to conform to ISO/IEC 9945-2:1993.  Please refer to the MacOS
man pages for details in regards to 'time'.


Model Identifier:	MacBookPro18,1

Chip:	Apple M1 Pro

Total Number of Cores:	10 (8 performance and 2 efficiency)

Memory:	32 GB

System Firmware Version:	7459.101.3

OS Loader Version:	7459.101.3

