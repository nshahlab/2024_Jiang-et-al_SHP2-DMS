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
		myfile = open(filename[:-4] + '_WPD distance.txt', 'w')
		structure = parser.get_structure(filename[:-4], filename)
		for model in structure:
			for chain in model:
                                dist425to459 = chain[425+1]["CA"] - chain[459+1]["CA"]
                                dist423to429 = chain[423+1]["CZ3"] - chain[429+1]["CG"]
                                dist424to428 = chain[424+1]["CG"] - chain[428+1]["C"]
                                myfile.write(str(dist425to459) + "\t" + str(dist423to429) + "\t" + str(dist424to428) + "\n")
		myfile.close()
