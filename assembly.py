# Hybrid bacterial genome assembly
# Author: Ciriac CHARLES
# Organization: ANSES Maisons Alfort and INRAe Nouzilly 
# Version: 2
# Date: 27/08/2022

import glob
import subprocess
import os
import sys
import shutil
import re
import argparse
import pathlib
import logging
from Bio import SeqIO


def run_with_env(env, cmd):
    cmd = cmd.replace("'", "'")
    command = f"bash -ic 'conda activate --stack {env} && {cmd} && conda deactivate'"
    process = subprocess.run(command, shell=True)
    if process.returncode != 0:
        LOGGER.error("Runtime error occured when running command:")
        LOGGER.error(cmd)
        sys.exit(process)


def simulated_reads():
    # Simulated reads
    os.makedirs(f"{PATH}/readsimulate{INPUT_DIR}", exist_ok=True)
    
    if not os.path.exists(f"{PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}shortread_R1.fq"):
        art_illumina_cmd = f"art_illumina \
            -i {PATH}/{INPUT_DIR}.fasta -f 50 -l 250 -p -na -m 300 -s 10 -ss MSv3 \
            -o {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}shortread_R"

        run_with_env("~/envs/art", art_illumina_cmd)

    if not os.path.exists(f"{PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz"):
        badread_cmd = f"badread simulate \
            --reference {PATH}/{INPUT_DIR}.fasta --quantity 50x --error_model random --qscore_model ideal \
            --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --identity 95,100,4 | \
            gzip > {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz"

        run_with_env("~/envs/Badread", badread_cmd)


def long_read_assemblies():
    # Long read assemblies (Unicycler, Raven et Flye)
    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/12longreads/assembly_12dir"):
        trycycler_cmd = f"trycycler subsample \
            -t 1 -r {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz \
            -o {PATH}/trycycler{INPUT_DIR}/12longreads/ --count 12 --genome_size 4366000"

        run_with_env("~/envs/trycycler", trycycler_cmd)

    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/12longreads/assembly_04dir"):
        for i in range(1, 5):
            n = "{:>02}".format(i)
            flye_cmd = f"flye \
                --nano-raw {PATH}/trycycler{INPUT_DIR}/12longreads/sample_{n}.fastq \
                --threads {THREADS} \
                --out-dir {PATH}/trycycler{INPUT_DIR}/12longreads/assembly_{n}dir"

            run_with_env("~/envs/flye", flye_cmd)

    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/12longreads/assembly_08.fasta"):
        for i in range(5, 9):
            n = "{:>02}".format(i)
            raven_cmd = f"raven \
                {PATH}/trycycler{INPUT_DIR}/12longreads/sample_{n}.fastq \
                > {PATH}/trycycler{INPUT_DIR}/12longreads/assembly_{n}.fasta"
            run_with_env("", raven_cmd)

    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/12longreads/assembly_12dir"):
        for i in range(9, 13):
            n = "{:>02}".format(i)
            raven_cmd = f"unicycler \
                -l {PATH}/trycycler{INPUT_DIR}/12longreads/sample_{n}.fastq \
                -o {PATH}/trycycler{INPUT_DIR}/12longreads/assembly_{n}dir "
            run_with_env("", raven_cmd)


def trycycler_consensus():
    # Consensus with trycycler
    os.makedirs(f"{PATH}/TRYCYCLER/", exist_ok=True)
    
    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/trycycler/cluster_001/"):
        # 	run_with_env("~/envs/trycycler", f"trycycler cluster --assemblies {PATH}/trycycler{INPUT_DIR}/12longreads/*.fasta {PATH}/trycycler{INPUT_DIR}/12longreads/*dir/*.fasta --reads {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz --out_dir {PATH}/trycycler{INPUT_DIR}/trycycler/")
        trycyler_cmd = f"trycycler cluster \
            --assemblies {PATH}/trycycler{INPUT_DIR}/12longreads/*.fasta {PATH}/trycycler{INPUT_DIR}/12longreads/*dir/*.fasta \
            --reads {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz \
            --out_dir {PATH}/trycycler{INPUT_DIR}/trycycler/"

        LOGGER.error("You need to complete the trycycler consensus by yourself. Go in the directory and run")
        LOGGER.error(trycyler_cmd)
        # 	sys.exit()

    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/trycycler/final_consensus.fasta"):
        if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/trycycler/cluster_002"):
            cluster = "cluster_001"
        else:
            cluster = input("Pick a cluster to use >")

        if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/trycycler/{cluster}/7_final_consensus.fasta"):
            trycycler_reconcile_cmd = f"trycycler reconcile \
                --linear --reads {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz \
                --cluster_dir {PATH}/trycycler{INPUT_DIR}/trycycler/{cluster}"
            run_with_env("~/envs/trycycler", trycycler_reconcile_cmd)

            trycycler_msa_cmd = f"trycycler msa --cluster_dir {PATH}/trycycler{INPUT_DIR}/trycycler/{cluster}"
            run_with_env("~/envs/trycycler", trycycler_msa_cmd)

            trycycler_partition_cmd = f"trycycler partition \
                --reads {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz \
                --cluster_dir {PATH}/trycycler{INPUT_DIR}/trycycler/{cluster}"
            run_with_env("~/envs/trycycler", trycycler_partition_cmd)

            trycycler_consensus_cmd = f"trycycler consensus --cluster_dir {PATH}/trycycler{INPUT_DIR}/trycycler/{cluster}"
            run_with_env("~/envs/trycycler", trycycler_consensus_cmd)

            shutil.copy(f"{PATH}/trycycler{INPUT_DIR}/trycycler/{cluster}/7_final_consensus.fasta", f"{PATH}/trycycler{INPUT_DIR}/final_consensus.fasta")


def medaka():
    # Medaka step
    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/medaka/consensus.fasta"):
        medaka_cmd = f"medaka_consensus \
            -i {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}longread.fastq.gz \
            -d {PATH}/trycycler{INPUT_DIR}/final_consensus.fasta \
            -o {PATH}/trycycler{INPUT_DIR}/medaka/ \
            -m r941_min_sup_g507 -t 8"

        run_with_env("~/envs/medaka", medaka_cmd)


def pilon():
    # Pilon step
    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/pilon_round4.fasta"):
        os.makedirs(f"{PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/", exist_ok=True)

        correction = 0
        current_assembly = f"{PATH}/trycycler{INPUT_DIR}/medaka/consensus.fasta"
        round = 1

        while correction < 2: 
            current_mapping = f"{PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/mapping_pilon{round}.sorted.bam"
            next_assembly = f"{PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/pilon_round{round+1}"
            run_with_env("~/envs/pilon", f"bwa index {current_assembly}")
            run_with_env("~/envs/pilon", f"bwa mem {current_assembly} {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}shortread_R1.fq {PATH}/readsimulate{INPUT_DIR}/{INPUT_DIR}shortread_R2.fq | samtools view - -Sb | samtools sort - -@14 -o {current_mapping}")
            run_with_env("~/envs/pilon", f"samtools index {current_mapping}")
            run_with_env("~/envs/pilon", f"java -jar pilon-1.24.jar  --genome {current_assembly} --fix all --changes --frags {current_mapping} --output {next_assembly} | tee {next_assembly}.pilon")
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

        shutil.copy(current_assembly, f"{PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/pilon_round.fasta")


def circlator():
    # Circularization with circlator
    os.makedirs(f"{PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/", exist_ok=True)

    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp.fasta"):
        run_with_env("", f"circlator fixstart --min_id 70 {PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/pilon_round.fasta {PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp")
        circlator_cmd = f"circlator fixstart \
            --min_id 70 \
            {PATH}/trycycler{INPUT_DIR}/pilon{INPUT_DIR}/pilon_round.fasta \
            {PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp"

        assembly = SeqIO.read(f"{PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp.fasta", "fasta")
        assembly.id = INPUT_DIR
        assembly.description = ""
        with open(f"{PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp_renamed.fasta", "w") as fasta_out:
            SeqIO.write(assembly, fasta_out, "fasta")


def assembly_analysis():
    # Assembly analysis
    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/quast/report.txt"):
        quast_cmd = f"quast.py \
            {PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp_renamed.fasta \
            -o {PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/quast"

        run_with_env("", quast_cmd)


def prokka():
    # Prokka with prodigal
    if not os.path.exists(f"{PATH}/prodigal.trn"):
        run_with_env("~/envs/prokka", f"prodigal -i {PATH}/Mb3601.gbk -t {PATH}/prodigal.trn")
    
    if not os.path.exists(f"{PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/prokka"):
        prokka_cmd = f"prokka \
            {PATH}/trycycler{INPUT_DIR}/circlatorstart/{INPUT_DIR}/{INPUT_DIR}exp_renamed.fasta \
            --prodigaltf prodigal.trn --protein {PATH}/Mb3601.gbk --prefix {INPUT_DIR} --force \
            --outdir {PATH}/circlatorstart/{INPUT_DIR}/"

        run_with_env("~/envs/prokka", prokka_cmd)


if __name__ == "__main__":
    # Argument parsers
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", help="Input directory", type=str, required=True)
    parser.add_argument("-p", "--path", help="Base path to be used (default: cwd {})".format(pathlib.Path().resolve()), type=str, default=pathlib.Path().resolve())
    parser.add_argument("-t", "--threads", help="Number of threads to use (default: 8)", type=int, default=8)
    parser.add_argument("--bypass-sr", help="Bypass simulated reads (default: False)", action='store_true', default=False)
    parser.add_argument("-v", "--verbosity", help="Output verbosity (default: 1)", type=int, choices=[0, 1, 2], default=1)
    cli_args = parser.parse_args()

    # Extracting GLOBALS
    THREADS = cli_args.threads
    PATH = cli_args.path
    VERBOSITY = {0: logging.ERROR, 1: logging.INFO, 2: logging.DEBUG}[cli_args.verbosity]
    INPUT_DIR = cli_args.input_dir

    # Setting up logger
    logging.basicConfig(level=VERBOSITY, format="[%(levelname)s] - %(message)s")
    LOGGER = logging.getLogger(__name__)
    LOGGER.debug(f"Log level set to {cli_args.verbosity}")

    # Running tools
    steps = [
        ["Simulated Reads", simulated_reads],
        ["Long read assemblied", long_read_assemblies],
        ["Trycycler consensus", trycycler_consensus],
        ["Medaka", medaka],
        ["Pilon", pilon],
        ["Circlator", circlator],
        ["Assembly analysis", assembly_analysis],
        ["Prokka", prokka],
    ]
    for step in steps:
        if cli_args.bypass_sr and step[1] == simulated_reads:
            continue

        LOGGER.info(f"Executing step: {step[0]}")
        step[1]()
        LOGGER.info(f"Completed step: {step[0]}")
