# ARACNe3

## ARACNe3 (Algorithm for the Reconstruction of Accurate Cellular Networks ver. 3) in C++

This is a prototype of the C++ implementation of ARACNe3 currently under development.  ARACNe3 presents multiple computational and theoretical improvements to the existing ARACNe-AP algorithm formulated in Java.  


Please contact Aaron T. Griffin and Lukas J. Vlahos for questions regarding this project.


Aaron T. Griffin - atg2142@cumc.columbia.edu 

Lukas J. Vlahos - lv2395@cumc.columbia.edu 

## Running the Most Recent Executable

This codebase is currently under development, but the most up-to-date ARACNe3 executable can still be run on any platform.  The executable must be compiled using source files in this directory, but only standard libraries are used, so this can be done by linking all object file dependencies with the ARACNe3 target executable, which can be run using the following commandline arguments: 


`./ARACNe3 /path/to/regulators.txt /path/to/gexpmatrix.txt`

The program will output an 'output.txt' that contains the regulator-target MI values estimated via Adaptive Partitioning (APMI), and the associated p-value determined from a null distribution of the APMI between 1 million shuffled gene-expression marginals.  Currently, this is computationally slow, can only run on one CPU core, and will take approximately 10 minutes on a computer with the specifications at the bottom of the page.   

## Current Developments:
 - Multithreading/non-multithreading option using standard library
 - **IN PROGRESS** Make null MI value independent of rest and stored on disk for ### subnet operations.  Maybe, check IF the file exists, if not, compute initNullMI
 - Remove any references on data types `<=4B`, as references instantiate pointer values which are at least 4B (and typically are 8B on 64-bit systems)
 - RETURN BY REFERENCE when applicable!!!
 - Optimize p-value calculation for each MI value
 - Low-level optimization and parallel for loop processing. Namely, minimizing heap allocation and using caches, as well as using the most efficient data structures required to store edge information (hashmaps, linked lists, adjacency matrices, etc.)
 
 Plan of action: Overhaul of Low-level algorithms -> Rewriting codebase with highest efficiency class design and matured codebase -> Implement multithreading 

## Tracking Progress

Whenever a significant change is made to an existing module of this program, test results are noted in the git commit notes.  The results reflect runtime performance of the most recent version of ARACNe3 or the given executable.  Tests were run on a computer with the specifications listed below.

## System Information
| System Feature | Value |
| :----: | :----: |
| Model Identifier: | MacBookPro18,1 |
| Chip: | Apple M1 Pro |
| Total Number of Cores: | 10 (8 performance and 2 efficiency) |
| Memory: | 32 GB |
| System Firmware Version: |7459.101.3 |
| OS Loader Version: | 7459.101.3 |
