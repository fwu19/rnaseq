#!/usr/bin/env bash

# module load R/4.1.2-foss-2020b
# export R_LIBS=/fh/fast/_SR/Genomics/user/fwu/R/x86_64-pc-linux-gnu-library/4.1 
module load Python/3.8.2-GCCcore-9.3.0
export PYTHONPATH=~/.local/lib/python3.8/site-packages/:${PATHONPATH}
export PATH=/home/fwu/.local/bin/:$PATH # pip3 install RSeQC --user

sampleId=$1; shift
rseqcBed=$1; shift
txBed=$1; shift
geneBed=$1; shift

for sampleBam in "$@"; do
    geneBody_coverage.py -i $sampleBam -r $rseqcBed -o $sampleId.GBC 2>$sampleId.log
    junction_annotation.py -i $sampleBam -o $sampleId.JA -r $txBed 2>>$sampleId.log
    read_duplication.py -i $sampleBam -o $sampleId.RD 2>>$sampleId.log
    read_distribution.py -i $sampleBam -r $txBed >$sampleId.RDis.read_distribution.txt 2>>$sampleId.log
    #inner_distance.py -i $sampleBam -o $sampleId.ID -r $geneBed 2>>$sampleId.log
done
	
