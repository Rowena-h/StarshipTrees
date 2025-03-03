#!/bin/bash
#SBATCH -p ei-short				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 10G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

out_dir=sourmash

mkdir -p ${out_dir}

realpath mashtree/starship_elements_big_split/* | awk 'BEGIN{FS=","; OFS=","} {print $1, $1}' | awk -F"," 'sub(".*/","",$1)' OFS="," | sed '1i name,genome_filename,protein_filename' > ${out_dir}/big_files.csv

#Big tree

singularity exec ~/programmes/sourmash/sourmash.img sourmash sketch fromfile \
	${out_dir}/big_files.csv \
	-p dna,k=21,scaled=10 \
	-o ${out_dir}/big_sketches.zip

singularity exec ~/programmes/sourmash/sourmash.img sourmash compare \
	${out_dir}/big_sketches.zip \
	--estimate-ani \
	--distance-matrix \
	--csv ${out_dir}/big_dist.csv

#Curated tree

realpath mashtree/starship_elements_conservative_split/* | awk 'BEGIN{FS=","; OFS=","} {print $1, $1}' | awk -F"," 'sub(".*/","",$1)' OFS="," | sed '1i name,genome_filename,protein_filename' > ${out_dir}/curated_files.csv

singularity exec ~/programmes/sourmash/sourmash.img sourmash sketch fromfile \
       ${out_dir}/curated_files.csv \
       -p dna,k=21,scaled=10 \
       -o ${out_dir}/curated_sketches.zip

singularity exec ~/programmes/sourmash/sourmash.img sourmash compare \
	${out_dir}/curated_sketches.zip \
	--estimate-ani \
	--distance-matrix \
	--csv ${out_dir}/curated_dist.csv
