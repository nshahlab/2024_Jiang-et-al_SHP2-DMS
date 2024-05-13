#!/usr/local/bin/python

import subprocess
import os
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PICIO, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.internal_coords import *
import string

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()
parser = PDBParser()
chainIDs = list(string.ascii_uppercase)

for filename in os.listdir(directory):
	if filename.endswith(".pdb"):
		myfile = open(filename[:-4] + '_R111-E249.txt', 'w')
		structure = parser.get_structure(filename[:-4], filename)
		for model in structure:
			for chain in model:
				dist01 = chain[111+1]["NH1"] - chain[249+1]["OE1"]
				dist02 = chain[111+1]["NH1"] - chain[249+1]["OE2"]
				dist03 = chain[111+1]["NH2"] - chain[249+1]["OE1"]
				dist04 = chain[111+1]["NH2"] - chain[249+1]["OE2"]
				distList = [dist01, dist02, dist03, dist04]
				myfile.write(str(dist01) + "\t" + str(dist02) + "\t" + str(dist03) + "\t" + str(dist04) + "\t" + str(min(distList)) + "\n")
		myfile.close()
