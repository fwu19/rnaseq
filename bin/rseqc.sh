#!/usr/bin/env bash

id=$1; shift
bam=$1; shift
bed=$1; shift

head -n 5000 ${bed} > subset.bed
log=${id}.log
geneBody_coverage.py -i $bam -r subset.bed -o ${id}.GBC 2>$log
junction_annotation.py -i $bam -r $bed -o ${id}.JA 2>>$log
read_duplication.py -i $bam -o ${id}.RD 2>>$log
read_distribution.py -i $bam -r $bed >${id}.RD.read_distribution.txt 2>>$log
#inner_distance.py -i $bam -o ${id}.ID -r $bed 2>>$log

	
