General notes:

These scripts were generated to process data for a manuscript entitled, "Deep mutational scanning of the multi-domain phosphatase SHP2 reveals mechanisms of regulation and pathogenicity" (https://doi.org/10.1101/2024.05.13.593907).
The script have been tested on Python 3.9.12 and run from the MacOS 12.7 terminal. Most of the scripts require BioPython and/or NumPy.
The scripts are divided into two subfolders, based on the type of data being analyzed: deep sequencing data from deep mutational scanning experiments and molecular dynamics (MD) trajectories.

Deep sequencing:

batch-flash.py - runs a program called FLASH for paired-end read merging (https://ccb.jhu.edu/software/FLASH/). This script will merge all of the paired reads in a directory, with R2 as the forward read and R1 as the reverse read.

batch-cutadapt.py - runs a program called CutAdapt (https://cutadapt.readthedocs.io/en/stable/) for adapter trimming after paired-end read merging using FLASH. For this script, you will need to change adapter sequence, depending on which SHP2 tile you are trimming data from. Sequences can be found in the file named "Cutadapt sequences.txt".

SHP2-AA-matrix.py - takes trimmed files that are the output from CutAdapt and counts the frequency of each mutation, at the amino acid level, within the specified tile. This script requires a tile sequence file "SHP2_AA_Tile_X.fa", where X should be replaced by the tile number you are analyzing. The input file name needs to be specified in the script.

MD analysis:

Each of the Python scripts in this file use the Biopython PDB module. These scripts will run in batch on all of the PDB files in a folder. If these files are MD trajectories, where each state is a different frame, the script will iterate over the states and calculate the desired parameter at each frame.
