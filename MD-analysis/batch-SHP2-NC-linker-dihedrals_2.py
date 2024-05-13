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
		myfile = open(filename[:-4] + '_NC-linker-dihedrals.txt', 'w')
		structure = parser.get_structure(filename[:-4], filename)
		for model in structure:
			for chain in model:
				phi103 = (calc_dihedral(chain[102+1]["C"].get_vector(),chain[103+1]["N"].get_vector(),chain[103+1]["CA"].get_vector(),chain[103+1]["C"].get_vector()))*180/np.pi
				psi103 = (calc_dihedral(chain[103+1]["N"].get_vector(),chain[103+1]["CA"].get_vector(),chain[103+1]["C"].get_vector(),chain[104+1]["N"].get_vector()))*180/np.pi
				phi104 = (calc_dihedral(chain[103+1]["C"].get_vector(),chain[104+1]["N"].get_vector(),chain[104+1]["CA"].get_vector(),chain[104+1]["C"].get_vector()))*180/np.pi
				psi104 = (calc_dihedral(chain[104+1]["N"].get_vector(),chain[104+1]["CA"].get_vector(),chain[104+1]["C"].get_vector(),chain[105+1]["N"].get_vector()))*180/np.pi
				phi105 = (calc_dihedral(chain[104+1]["C"].get_vector(),chain[105+1]["N"].get_vector(),chain[105+1]["CA"].get_vector(),chain[105+1]["C"].get_vector()))*180/np.pi
				psi105 = (calc_dihedral(chain[105+1]["N"].get_vector(),chain[105+1]["CA"].get_vector(),chain[105+1]["C"].get_vector(),chain[106+1]["N"].get_vector()))*180/np.pi
				phi106 = (calc_dihedral(chain[105+1]["C"].get_vector(),chain[106+1]["N"].get_vector(),chain[106+1]["CA"].get_vector(),chain[106+1]["C"].get_vector()))*180/np.pi
				psi106 = (calc_dihedral(chain[106+1]["N"].get_vector(),chain[106+1]["CA"].get_vector(),chain[106+1]["C"].get_vector(),chain[107+1]["N"].get_vector()))*180/np.pi
				myfile.write(str(phi103) + "\t" + str(psi103) + "\t" + str(phi104) + "\t" + str(psi104) + "\t" + str(phi105) + "\t" + str(psi105) + "\t" + str(phi106) + "\t" + str(psi106) + "\n")
		myfile.close()
