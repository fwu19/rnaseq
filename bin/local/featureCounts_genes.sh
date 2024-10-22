#!/usr/bin/env bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=fwu@fredhutch.org
#SBATCH --output=/fh/scratch/delete30/_SR/Genomics/fwu/%A-%a.%x.out
#SBATCH --error=/fh/scratch/delete30/_SR/Genomics/fwu/%A-%a.%x.err
#SBATCH --nodes=1
#SBATCH --mem 36G
#SBATCH -n 6

###################################################################################
#### Count reads for genes using featureCounts
###################################################################################
if [[ $# < 1 ]]; then
	echo "This script count reads covering gene models using featureCounts."
	echo -e "USAGE: bash `basename $0` <params> <sample IDs separated by space>"
	exit 1;
fi

params=$1; shift
source $params
[[ -z $analysisDir ]] && { echo '$analysisDir is missing'; exit 1; }
[[ -z $gtfFile ]] && { echo '$gtfFile is missing'; exit 1; }
[[ -z $strand ]] && { echo '$strand is missing'; exit 1; }

## Read samples
samples="$1"; shift # a file to store sample IDs, one per line
if [[ -f $samples ]]; then
	IFS=$'\n\r' read -d '' -r -a samplelist < $samples # read csv file
else
	IFS=' ' read -r -a samplelist <<< "$samples" # read the string
fi	

# load modules and conda environment
source ~/.bashrc
module load Subread/2.0.0-GCC-8.3.0

## cd analysis directory
workDir=$analysisDir/featureCounts
[[ -d $workDir ]] || mkdir -p $workDir
cd $workDir

## Process one sample
IFS=',' read -r -a arr <<< "${samplelist[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"        
sample=${arr[0]}

# count fragments for exons only for paired-end reads, do not allow multi-overlapping, without -O --fraction 
## for PE, use -B from 05/13/2020 and stop later
if [[ $readType == PE ]]; then
	featureCounts -T ${core:-1} -p -C --minOverlap 1 -a ${gtfFile4featureCounts:-$gtfFile} -F GTF -t ${featureType:-exon} -g ${attrType:-gene_id} -s $strand -o $sample.exonicFragmentCounts.txt $analysisDir/alignments/$sample.bam 2>$sample.err.log; 
else
        featureCounts -T ${core:-1} --minOverlap 1 -a ${gtfFile4featureCounts:-$gtfFile} -F GTF -t ${featureType:-exon} -g ${attrType:-gene_id} -s $strand -o $sample.exonicReadCounts.txt $analysisDir/alignments/$sample.bam 2>$sample.err.log; 
fi

