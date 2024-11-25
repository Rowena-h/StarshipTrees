#!/bin/bash
#SBATCH -p ei-medium				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 10					# number of cores
#SBATCH --mem 100G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

#Extract elements into separate fasta files

#seqkit 0.12.0
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6

seqkit split \
	--by-id Gluck-Thaler2024/sequences/sequences/mycodb.final.starships.fna \
	-O mashtree/starship_elements_big_split

rename .fa .fasta mashtree/starship_elements_big_split/*
rename .fna .fasta mashtree/starship_elements_big_split/*

cp mashtree/starship_elements_split/gaeumannomyces* mashtree/starship_elements_big_split

#Run mashtree
singularity exec ~/programmes/mashtree/mashtree.img mashtree_bootstrap.pl \
	--reps 100 \
	--numcpus ${SLURM_CPUS_PER_TASK} \
	 mashtree/starship_elements_big_split/* \
	-- > mashtree/starship_tree_big.bootstrap.tre

#Repeat only with elements with flanking direct repeats
cp -r mashtree/starship_elements_big_split mashtree/starship_elements_big_conservative_split

cd mashtree/starship_elements_big_conservative_split

for file in $(ls)
do
	new=$(echo $file | sed 's/__.*/\.fasta/')
	mv $file $new
done

comm -2 -3 <(ls) <(sort ../../flanked) | xargs rm

cp mashtree/starship_elements_conservative_split/gaeumannomyces.elements* mashtree/starship_elements_big_conservative_split

singularity exec ~/programmes/mashtree/mashtree.img mashtree_bootstrap.pl \
        --reps 100 \
        --numcpus ${SLURM_CPUS_PER_TASK} \
         mashtree/starship_elements_big_conservative_split/* \
        -- > mashtree/starship_tree_big_conservative.bootstrap.tre
