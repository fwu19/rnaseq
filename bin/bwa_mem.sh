#!/usr/bin/env bash

cpus=$1;shift
bwa_index=$1; shift
read1=$1; shift
read2=$1; shift
out_prefix=$1; shift

bwa mem -M -t ${cpus} ${bwa_index}/$(ls ${bwa_index} | egrep "bwt$" | sed "s/.bwt$//") <(zcat $read1) <(zcat $read2) | samtools sort -@ ${cpus} - >${out_prefix}.bam 
samtools index ${out_prefix}.bam