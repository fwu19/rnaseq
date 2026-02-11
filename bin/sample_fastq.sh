#!/usr/bin/env bash

fq1=$1; shift
fq2=$1; shift
out_prefix=$1; shift
n_reads=$1; shift
let n_lines=$n_reads*4

(zcat ${fq1} | head -n ${n_lines} | gzip -c > ${out_prefix}_R1.head_${n_reads}.fastq.gz); 

(zcat ${fq2} | head -n ${n_lines} | gzip -c > ${out_prefix}_R2.head_${n_reads}.fastq.gz); 
status=$?
[ $status -eq 141 ] && status=0
exit $status
