#!/usr/bin/env bash

sample_id=$1; shift

cat $(find ./ | egrep "S[0-9]+(_L[0-9]+)?_R1_" | sort| tr '\n' ' ' ) > ${sample_id}_merged_R1.fastq.gz
cat $(find ./ | egrep "S[0-9]+(_L[0-9]+)?_R2_" | sort| tr '\n' ' ' ) > ${sample_id}_merged_R2.fastq.gz
