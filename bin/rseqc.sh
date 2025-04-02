#!/usr/bin/env bash

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
	
