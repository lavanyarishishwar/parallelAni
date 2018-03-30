# parallelAni
Embarrassingly parallel implementation of ANI from whole genomes


Parallel script for computing pairwise ANIm values.
Requires mummer to be installed.

Usage instructions:
getANIm-parallel.pl [-out <output ANI matrix. Default:ani.txt>] 
                    [-ext <extension of the fasta files. Default:fasta>] 
                    [-threads <number of threads. Default: 10>] 
                    [-folder <Folder in which to look for file. Default: ./>] 
                    [-distance <FLAG. Prints distances (100-ani) instead of similarity>] 
                    [-help <FLAG. Prints the help>]

