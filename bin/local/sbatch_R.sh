#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.out
#SBATCH --error=/hpc/temp/_SR/Genomics/users/fwu/slurm/%A-%a.%x.err

module purge all
module load fhR/4.1.2-foss-2021b
export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 # 

Rscript $@