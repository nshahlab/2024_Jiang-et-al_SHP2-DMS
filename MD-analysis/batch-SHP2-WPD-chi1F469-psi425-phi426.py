#!/usr/local/bin/python

import subprocess
import os
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PICIO, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.internal_coords import *
from Bio.PDB.vectors import calc_dihedral
import string

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()
parser = PDBParser()
chainIDs = list(string.ascii_uppercase)

for filename in os.listdir(directory):
	if filename.endswith(".pdb"):
		myfile = open(filename[:-4] + '_WPD angles.txt', 'w')
		structure = parser.get_structure(filename[:-4], filename)
		for model in structure:
			for chain in model:
				chi1F469 = (calc_dihedral(chain[469+1]["N"].get_vector(),chain[469+1]["CA"].get_vector(),chain[469+1]["CB"].get_vector(),chain[469+1]["CG"].get_vector()))*180/np.pi
				psi425 = (calc_dihedral(chain[425+1]["N"].get_vector(),chain[425+1]["CA"].get_vector(),chain[425+1]["C"].get_vector(),chain[426+1]["N"].get_vector()))*180/np.pi
				phi426 = (calc_dihedral(chain[425+1]["C"].get_vector(),chain[426+1]["N"].get_vector(),chain[426+1]["CA"].get_vector(),chain[426+1]["C"].get_vector()))*180/np.pi
				
				myfile.write(str(chi1F469) + "\t" + str(psi425) + "\t" + str(phi426) + "\n")
		myfile.close()
