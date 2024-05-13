#!/usr/local/bin/python

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
from collections import OrderedDict

# Import data

wtSeq = SeqIO.read(open("SHP2_AA_Tile_1.fa"), "fasta")
aaList = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"]
wtSeqList = list(wtSeq.seq)



# Enumerate variants

variants = {str(wtSeq.seq):0}
for position, residue in enumerate(wtSeqList):
	remaining = list(set(aaList) - set(residue))
	for aa in remaining:
		variant = ''.join(wtSeqList[0:position]) + ''.join(aa) + ''.join(wtSeqList[position+1:])
		variants[variant] = 0

# Translate nucleotide sequence and count variants

for nuc in SeqIO.parse("D706-507_S47_L001.trimmed.fastq", "fastq"):
	protein = SeqRecord(seq = nuc.seq.translate(to_stop=False), id = "trans_" + nuc.id, description = "translation of variant")
	if str(protein.seq) in variants:
		variants[str(protein.seq)] += 1

# Generate count matrix

countMatrix = numpy.zeros(shape=(len(aaList),len(wtSeqList)))

for pos in range(len(wtSeqList)):
	for aa in range(len(aaList)):
		sequence = ''.join(wtSeqList[0:pos]) + ''.join(aaList[aa]) + ''.join(wtSeqList[pos+1:])
		countMatrix[aa][pos] = variants[sequence]

numpy.savetxt("D706-507_S47_L001.count_matrix.txt", countMatrix, fmt = '%i', delimiter='\t')
