#!/usr/bin/env bash

sample_id=$1; shift
indir=$1; shift
cat $(find $indir/ | egrep "S[0-9]+(_L[0-9]+)?_R1_" | tr '\n' ' ' ) > ${sample_id}_merged_R1.fastq.gz
cat $(find $indir/ | egrep "S[0-9]+(_L[0-9]+)?_R2_" | tr '\n' ' ' ) > ${sample_id}_merged_R2.fastq.gz
