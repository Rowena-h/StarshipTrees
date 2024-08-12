#!/bin/bash
#SBATCH -p ei-short				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 1G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

strain=$(awk '{print $1}' strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
proteins=../data/ncbi_data/proteins/gaeumannomyces/${strain}_EIv1.release.gff3.pep.repisoform.fasta

source blast-2.10

makeblastdb \
	-in ${proteins} \
	-dbtype prot

for prot in spok2 toxa NLR PLP FRE DUF3723 phi-base_240801
do

	blastp 	-query blasts/${prot}.fa \
		-db ${proteins} \
		-outfmt "6 qseqid sseqid evalue bitscore pident length" \
		-evalue 1e-25 \
		-out blasts/${strain}_${prot}.tsv \
		-num_threads ${SLURM_CPUS_PER_TASK}

done

for nuc in spok1 spok3 spok4
do

	blastx 	-query ${nuc}.fa \
		-db ${proteins} \
		-outfmt "6 qseqid sseqid evalue bitscore pident length" \
		-evalue 1e-25 \
		-out blasts/${strain}_${nuc}.tsv \
		-num_threads ${SLURM_CPUS_PER_TASK}

done
