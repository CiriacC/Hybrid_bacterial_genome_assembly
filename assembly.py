#Ciriac CHARLES INRAe Nouzilly 2021


import glob
import subprocess
import os
import sys
import shutil
import re
from Bio import SeqIO


def run_with_env(env, cmd):
    cmd= cmd.replace("'", "\'")
    command= f"bash -ic 'conda activate --stack {env} && {cmd} && conda deactivate'"
    process = subprocess.run(command, shell=True)
    if process.returncode != 0:
        sys.exit(process)

if __name__ == "__main__":
    THREADS=8
    PATH="/home/ccharles/Documents"
    input_dir= sys.argv[1]


#Simulated reads

if not os.path.exists(f"{PATH}//readsimulate{input_dir}"):
	os.makedirs(f"{PATH}/readsimulate{input_dir}", exist_ok=True)
if not os.path.exists(f"{PATH}/readsimulate{input_dir}/{input_dir}shortread_R1.fq"):
	run_with_env("~/envs/art", f"art_illumina -i {PATH}/{input_dir}.fasta -f 50 -l 250 -p -na -m 300 -s 10 -ss MSv3 -o {PATH}/readsimulate{input_dir}/{input_dir}shortread_R")
if not os.path.exists(f"{PATH}/readsimulate{input_dir}/{input_dir}longread.fastq.gz"):
	run_with_env("~/envs/Badread", f"badread simulate --reference {PATH}/{input_dir}.fasta --quantity 50x --error_model random --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --identity 95,100,4 | gzip > {PATH}/readsimulate{input_dir}/{input_dir}longread.fastq.gz")

#Long read assemblies (Unicycler, Raven et Flye)
threads=16  # change as appropriate for your system

if not os.path.exists (f"{PATH}/trycycler{input_dir}/12longreads/assembly_12dir"):
    run_with_env("~/envs/trycycler", f"trycycler subsample -t 1 -r {PATH}/MinION/trim/{input_dir}_trimmed.fastq.gz -o {PATH}/trycycler{input_dir}/12longreads/ --count 12 --genome_size 4366000")

if not os.path.exists(f"{PATH}/trycycler{input_dir}/12longreads/assembly_04dir"):
	current_flye = (f"flye --nano-raw {PATH}/trycycler{input_dir}/12longreads/")
	current_flye2 = (f"--threads 8 --out-dir {PATH}/trycycler{input_dir}/12longreads/")
	run_with_env("~/envs/flye", f"{current_flye}sample_01.fastq {current_flye2}assembly_01dir")
	run_with_env("~/envs/flye", f"{current_flye}sample_02.fastq {current_flye2}assembly_02dir")
	run_with_env("~/envs/flye", f"{current_flye}sample_03.fastq {current_flye2}assembly_03dir")
	run_with_env("~/envs/flye", f"{current_flye}sample_04.fastq {current_flye2}assembly_04dir")

if not os.path.exists(f"{PATH}/trycycler{input_dir}/12longreads/assembly_08.fasta"):
	current_raven = (f"raven {PATH}/trycycler{input_dir}/12longreads/")
	current_raven2 = (f"> {PATH}/trycycler{input_dir}/12longreads/")
	run_with_env("", f"{current_raven}sample_05.fastq {current_raven2}assembly_05.fasta")
	run_with_env("", f"{current_raven}sample_06.fastq {current_raven2}assembly_06.fasta")
	run_with_env("", f"{current_raven}sample_07.fastq {current_raven2}assembly_07.fasta")
	run_with_env("", f"{current_raven}sample_08.fastq {current_raven2}assembly_08.fasta")

if not os.path.exists(f"{PATH}/trycycler{input_dir}/12longreads/assembly_12dir"):
	current_unicycler = (f"unicycler -l {PATH}/trycycler{input_dir}/12longreads/")
	current_unicycler2 = (f"-o {PATH}/trycycler{input_dir}/12longreads/")
	run_with_env("~/envs/unicycler", f"{current_unicycler}sample_09.fastq {current_unicycler2}assembly_09dir")
	run_with_env("~/envs/unicycler", f"{current_unicycler}sample_10.fastq {current_unicycler2}assembly_10dir")
	run_with_env("~/envs/unicycler", f"{current_unicycler}sample_11.fastq {current_unicycler2}assembly_11dir")
	run_with_env("~/envs/unicycler", f"{current_unicycler}sample_12.fastq {current_unicycler2}assembly_12dir")

#consensus with trycycler

if not os.path.exists(f"{PATH}/trycycler{input_dir}/"):
	os.makedirs(f"{PATH}/TRYCYCLER/", exist_ok=True)
if not os.path.exists(f"{PATH}/trycycler{input_dir}/trycycler/cluster_001/"):
#	run_with_env("~/envs/trycycler", f"trycycler cluster --assemblies {PATH}/trycycler{input_dir}/12longreads/*.fasta {PATH}/trycycler{input_dir}/12longreads/*dir/*.fasta --reads {PATH}/MinION/trim/{input_dir}_trimmed.fastq.gz --out_dir {PATH}/trycycler{input_dir}/trycycler/")
	print("il faut le lancer à la main. Allez dans la directory pour le lancer")
	print(f"trycycler cluster --assemblies {PATH}/trycycler{input_dir}/12longreads/*.fasta {PATH}/trycycler{input_dir}/12longreads/*dir/*.fasta --reads {PATH}/MinION/trim/{input_dir}_trimmed.fastq.gz --out_dir {PATH}/trycycler{input_dir}/trycycler/")
#	sys.exit()

if not os.path.exists(f"{PATH}/trycycler{input_dir}/trycycler/final_consensus.fasta"):
	if not os.path.exists(f"{PATH}/trycycler{input_dir}/trycycler/cluster_002"):
		cluster = "cluster_001"
	else:
		cluster = input("quel cluster choisir?")

	if not os.path.exists(f"{PATH}/trycycler{input_dir}/trycycler/{cluster}/7_final_consensus.fasta"):
		run_with_env("~/envs/trycycler", f"trycycler reconcile --linear --reads {PATH}/MinION/trim/{input_dir}_trimmed.fastq.gz --cluster_dir {PATH}/trycycler{input_dir}/trycycler/{cluster}")
		run_with_env("~/envs/trycycler", f"trycycler msa --cluster_dir {PATH}/trycycler{input_dir}/trycycler/{cluster}")
		run_with_env("~/envs/trycycler", f"trycycler partition --reads {PATH}/MinION/trim/{input_dir}_trimmed.fastq.gz  --cluster_dir {PATH}/trycycler{input_dir}/trycycler/{cluster}")
		run_with_env("~/envs/trycycler", f"trycycler consensus --cluster_dir {PATH}/trycycler{input_dir}/trycycler/{cluster}")
		shutil.copy(f"{PATH}/trycycler{input_dir}/trycycler/{cluster}/7_final_consensus.fasta", f"{PATH}/trycycler{input_dir}/final_consensus.fasta")

#Medaka step

if not os.path.exists(f"{PATH}/trycycler{input_dir}/medaka/consensus.fasta"):
	run_with_env("~/envs/medaka", f"medaka_consensus -i {PATH}/MinION/trim/{input_dir}_trimmed.fastq.gz -d {PATH}/trycycler{input_dir}/final_consensus.fasta -o {PATH}/trycycler{input_dir}/medaka/ -m r941_min_sup_g507 -t 8")

##Pilon step

if not os.path.exists(f"{PATH}/trycycler{input_dir}/pilon{input_dir}/pilon_round4.fasta"):
	os.makedirs(f"{PATH}/trycycler{input_dir}/pilon{input_dir}/", exist_ok=True)

correction = 0
current_assembly = f"{PATH}/trycycler{input_dir}/medaka/consensus.fasta"
round = 1

while correction < 2: 
	current_mapping = f"{PATH}/trycycler{input_dir}/pilon{input_dir}/mapping_pilon{round}.sorted.bam"
	next_assembly = f"{PATH}/trycycler{input_dir}/pilon{input_dir}/pilon_round{round+1}"
	run_with_env("~/envs/pilon", f"bwa index {current_assembly}")
	run_with_env("~/envs/pilon", f"bwa mem {current_assembly} {PATH}/illumina_wgs_trim_{input_dir}_R1_001.fastq.gz {PATH}/illumina_wgs_trim_{input_dir}_R2_001.fastq.gz | samtools view - -Sb | samtools sort - -@14 -o {current_mapping}")
	run_with_env("~/envs/pilon", f"samtools index {current_mapping}")
	run_with_env("~/envs/pilon", f"pilon --genome {current_assembly} --fix all --changes --frags {current_mapping} --output {next_assembly} | tee {next_assembly}.pilon")
	round += 1
	current_assembly = next_assembly+".fasta"
	with open (f"{next_assembly}.pilon") as pilon_log: 
		data = pilon_log.read()
		match = re.search("Corrected ([0-9]+) snps; 0 ambiguous bases; corrected ([0-9]+) small insertions totaling 0 bases, ([0-9]+) small deletions totaling 0 bases",data) 
		if match: 
			snps, insertion, deletion = [int(i) for i in match.groups()]
			if snps or insertion or deletion: 
				continue
			else: 
				correction += 1

shutil.copy(current_assembly, f"{PATH}/trycycler{input_dir}/pilon{input_dir}/pilon_round.fasta")

#circularization with circlator

if not os.path.exists(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/"):
	os.makedirs(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/", exist_ok=True)
if not os.path.exists(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/{input_dir}exp.fasta"):
	run_with_env("~/envs/circlator", f"circlator fixstart --min_id 70 {PATH}/trycycler{input_dir}/pilon{input_dir}/pilon_round.fasta {PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/{input_dir}exp")

	assembly = SeqIO.read(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/{input_dir}exp.fasta", "fasta")
	assembly.id = input_dir
	assembly.description= ""
	with open(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/{input_dir}exp_renamed.fasta", "w") as fasta_out:
		SeqIO.write(assembly, fasta_out, "fasta")


#Assembly analysis

if not os.path.exists(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/quast/report.txt"):
	run_with_env("", f"quast.py {PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/{input_dir}exp_renamed.fasta -o {PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/quast")

#Prokka with prodigal

if not os.path.exists(f"{PATH}/prodigal.trn"):
	run_with_env("~/envs/prokka", f"prodigal -i {PATH}/Mb3601.gbk -t {PATH}/prodigal.trn")
if not os.path.exists(f"{PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/prokka"):
	run_with_env("~/envs/prokka", f"prokka {PATH}/trycycler{input_dir}/circlatorstart/{input_dir}/{input_dir}exp_renamed.fasta --prodigaltf prodigal.trn --protein {PATH}/Mb3601.gbk --prefix {input_dir} --force --outdir {PATH}/circlatorstart/{input_dir}/ ")
