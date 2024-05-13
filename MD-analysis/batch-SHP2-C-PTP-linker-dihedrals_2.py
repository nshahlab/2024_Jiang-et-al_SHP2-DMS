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
		myfile = open(filename[:-4] + '_C-PTP-linker-dihedrals.txt', 'w')
		structure = parser.get_structure(filename[:-4], filename)
		for model in structure:
			for chain in model:
				phi216 = (calc_dihedral(chain[215+1]["C"].get_vector(),chain[216+1]["N"].get_vector(),chain[216+1]["CA"].get_vector(),chain[216+1]["C"].get_vector()))*180/np.pi
				psi216 = (calc_dihedral(chain[216+1]["N"].get_vector(),chain[216+1]["CA"].get_vector(),chain[216+1]["C"].get_vector(),chain[217+1]["N"].get_vector()))*180/np.pi
				phi217 = (calc_dihedral(chain[216+1]["C"].get_vector(),chain[217+1]["N"].get_vector(),chain[217+1]["CA"].get_vector(),chain[217+1]["C"].get_vector()))*180/np.pi
				psi217 = (calc_dihedral(chain[217+1]["N"].get_vector(),chain[217+1]["CA"].get_vector(),chain[217+1]["C"].get_vector(),chain[218+1]["N"].get_vector()))*180/np.pi
				phi218 = (calc_dihedral(chain[217+1]["C"].get_vector(),chain[218+1]["N"].get_vector(),chain[218+1]["CA"].get_vector(),chain[218+1]["C"].get_vector()))*180/np.pi
				psi218 = (calc_dihedral(chain[218+1]["N"].get_vector(),chain[218+1]["CA"].get_vector(),chain[218+1]["C"].get_vector(),chain[219+1]["N"].get_vector()))*180/np.pi
				phi219 = (calc_dihedral(chain[218+1]["C"].get_vector(),chain[219+1]["N"].get_vector(),chain[219+1]["CA"].get_vector(),chain[219+1]["C"].get_vector()))*180/np.pi
				psi219 = (calc_dihedral(chain[219+1]["N"].get_vector(),chain[219+1]["CA"].get_vector(),chain[219+1]["C"].get_vector(),chain[220+1]["N"].get_vector()))*180/np.pi
				myfile.write(str(phi216) + "\t" + str(psi216) + "\t" + str(phi217) + "\t" + str(psi217) + "\t" + str(phi218) + "\t" + str(psi218) + "\t" + str(phi219) + "\t" + str(psi219) + "\n")
		myfile.close()
