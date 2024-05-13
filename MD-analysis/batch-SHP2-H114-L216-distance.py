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
		myfile = open(filename[:-4] + '_H114-L216.txt', 'w')
		structure = parser.get_structure(filename[:-4], filename)
		for model in structure:
			for chain in model:
				dist01 = chain[114+1]["NE2"] - chain[216+1]["O"]
				myfile.write(str(dist01) + "\n")
		myfile.close()
