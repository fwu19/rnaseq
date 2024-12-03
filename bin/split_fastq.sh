#!/usr/bin/env bash

id=$1; shift
fq1=$1; shift
fq2=$1; shift
args="$1"; shift

seqkit split2 -1 $fq1 -2 $fq2 $args
r1=$(echo $fq1 | sed "s/.fastq.gz//" ).part_
r2=$(echo $fq2 | sed "s/.fastq.gz//" ).part_

: << 'NOT_RUN'
r1=$(echo $fq1 | sed "s/.fastq.gz//" ).part_
zcat $fq1 | split -d -l 100000000 - $r1
for i in $(find . | egrep $r1 ); do
    mv $i ${i}.fastq
    gzip ${i}.fastq
done

r2=$(echo $fq2 | sed "s/.fastq.gz//" ).part_
zcat $fq2 | split -d -l 100000000 - $r2
for i in $(find . | egrep $r2 ); do
    mv $i ${i}.fastq
    gzip ${i}.fastq
done
NOT_RUN


paste -d "," <(find ./ | egrep $r1 | sort | sed "s/.*part_/part_/g; s/.fastq.*$//" ) <(find $(realpath ./) | egrep $r1 | sort )  <(find $(realpath ./) | egrep $r2 | sort ) | \
awk -v id="${id}" 'BEGIN {OFS=FS=","} {print id,id"."$1,$2,$3}' >${id}.csv

