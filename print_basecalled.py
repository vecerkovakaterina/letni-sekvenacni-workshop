import pod5
import subprocess
import os

# TODO path to the first pod5 file

path_to_dorado_models="/home/lab33/researchers_night/dna_r10.4.1_e8.2_400bps_sup@v4.3.0/"
first_pod5 = "/var/lib/minknowOld/20241023_Artemia/Artemia/20241023_1453_MN35406_FAY54811_f133089e/pod5/FAY54811_f133089e_b6e97c76_3.pod5"
##first_pod5 = "/var/lib/minknow/data/20240927_S_pastoriani/S_pastoriani/20240927_1613_MN35406_FAY54794_83b7ce44/pod5/FAY54794_83b7ce44_e0878d97_0.pod5"

record_counter = 0

try:
	os.remove("tmp.pod5")
except OSError:
	pass

while True:
	with pod5.Reader(first_pod5) as reader:
		reads = list(reader.reads())
		records = len(reads)
		# if new record added to pod5 file
		if records > record_counter:
			record_counter += 1
			# write the newly added pod5 record to separate file
			with pod5.Writer("tmp.pod5") as writer:
				last_read = reads[record_counter].to_read()
				writer.add_read(last_read)
			# basecall
			##print("Basecalling the next read / Louskám další čtení...")
			command = f"dorado basecaller --min-qscore 25 {path_to_dorado_models} tmp.pod5 > tmp.bam"
			result = subprocess.run(command,
						shell=True,
						stdout=subprocess.PIPE,
						stderr=subprocess.PIPE,
						text=True)
			if result.returncode != 0:
				raise Exception(f"Failed to basecall {last_read.read_id}!")
			# bam to fasta
			command = "samtools fasta tmp.bam > tmp.fasta"
			result = subprocess.run(command,
						shell=True,
						stdout=subprocess.PIPE,
						stderr=subprocess.PIPE,
						text=True)
			if result.returncode != 0:
				raise Exception(f"Failed to convert to fasta {last_read.read_id}!")
			# print fasta file to console
			with open("tmp.fasta") as fasta_file:
				content = fasta_file.read()
				if len(content) > 0:
					print(content)
	os.remove("tmp.pod5")
	os.remove("tmp.bam")
	os.remove("tmp.fasta")
