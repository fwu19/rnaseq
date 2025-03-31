#!/usr/bin/env bash

sample_id=$1; shift
read1_dir=$1; shift
read2_dir=$1; shift

cat $(find ${read1_dir}/ | sort| tr '\n' ' ' ) > ${sample_id}_merged_R1.fastq.gz
cat $(find ${read2_dir}/ | sort| tr '\n' ' ' ) > ${sample_id}_merged_R2.fastq.gz
