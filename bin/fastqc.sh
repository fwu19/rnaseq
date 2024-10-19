#!/usr/bin/env bash
	
sample_id=$1; shift
indir=$1; shift
fq1=$(find $indir/ | egrep "merged_R1.fastq.gz|S[0-9]+(_L[0-9]+)?_R1_.*fastq.gz" )
fq2=$(find $indir/ | egrep "merged_R2.fastq.gz|S[0-9]+(_L[0-9]+)?_R2_.*fastq.gz" )
fastqc -o ./ --casava $fq1 $fq2 2>&1 | tee ${sample_id}.fastqc.log

