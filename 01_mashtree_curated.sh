#!/bin/bash
#SBATCH -p ei-medium				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 5					# number of cores
#SBATCH --mem 1G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

#Extract elements into separate fasta files

#seqkit 0.12.0
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6

seqkit split \
	--by-id ../GaeumannomycesGenomics/07_comparative_genomics/starfish/elementFinder/gaeumannomyces.elements.fna \
	-O mashtree/starship_elements_split

seqkit split \
	--by-id Gluck-Thaler2022/starships/Starships.fa \
	-O mashtree/starship_elements_split

rename .fa .fasta mashtree/starship_elements_split/*
rename .fna .fasta mashtree/starship_elements_split/*

#Run mashtree
singularity exec ~/programmes/mashtree/mashtree.img mashtree_bootstrap.pl \
	--reps 1000 \
	--numcpus ${SLURM_CPUS_PER_TASK} \
	mashtree/starship_elements_split/* \
	-- --min-depth 0 > mashtree/starship_tree.bootstrap.tre