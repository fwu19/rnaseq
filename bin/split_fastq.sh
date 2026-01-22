#!/usr/bin/env bash

id=$1; shift
fq1=$1; shift
fq2=$1; shift
size=$1; shift
args="$1"; shift

seqkit split2 -1 $fq1 -2 $fq2 -s $size $args
r1=$(echo $fq1 | sed "s/.fastq.gz//" ).part_
r2=$(echo $fq2 | sed "s/.fastq.gz//" ).part_

paste -d "," <(find ./ | egrep $r1 | sort | sed "s/.*part_/part_/g; s/.fastq.*$//" ) <(find $(realpath ./) | egrep $r1 | sort )  <(find $(realpath ./) | egrep $r2 | sort ) | \
awk -v id="${id}" 'BEGIN {OFS=FS=","} {print id,id"."$1,$2,$3}' >${id}.csv

