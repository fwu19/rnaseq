#!/usr/bin/env bash

sample_id=$1; shift
read1_dir=$1; shift
read2_dir=$1; shift

cat $(ls ${read1_dir}/*/* | sort| tr '\n' ' ' ) > ${sample_id}_merged_R1.fastq.gz
cat $(ls ${read2_dir}/*/* | sort| tr '\n' ' ' ) > ${sample_id}_merged_R2.fastq.gz
