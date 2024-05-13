import subprocess
import os

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()

for filename in os.listdir(directory):
	if filename.endswith(".extendedFrags.fastq"):
		trimmed = filename[0:-20] + ".trimmed.fastq"
		subprocess.run(['/opt/miniconda3/bin/cutadapt', '-g', 'NNNTGTTTAACTTTAAGAAGGAGATATACC...GATTTTACTCTCTCCGTTAGAAGNNN', '--discard-untrimmed', '-o', trimmed, filename])

