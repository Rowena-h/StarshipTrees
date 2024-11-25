#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 2                                    # number of cores
#SBATCH --mem 1GB                               # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

dir=raxmlng
aln=${dir}/captains_aln_trim.fasta
model=JTT+I+G4

#Convert file for RAxML-NG
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--parse \
	--msa ${aln} \
	--model ${model} \
	--prefix ${dir}/captains

#Run ML tree search and bootstrapping for <=1000 iterations
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--all \
	--msa ${aln} \
	--model ${model} \
	--prefix ${dir}/captains \
	--seed 2 \
	--threads ${SLURM_CPUS_PER_TASK} \
	--bs-trees autoMRE{1000}

#Check convergence
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--bsconverge \
	--bs-trees ${dir}/captains.raxml.bootstraps \
	--prefix ${dir}/captains_convergence_test \
	--seed 2
