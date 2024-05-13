import subprocess
import os

directory = subprocess.check_output(['pwd'], universal_newlines=True).rstrip()

for filename in os.listdir(directory):
	if filename.endswith("_R2_001.fastq"):
		filename2 = filename[0:-13] + "_R1_001.fastq"
		subprocess.run(['/Users/neelshah/Applications/FLASH2-2.2.00/flash2', filename, filename2, '-o', filename[0:-13]])
