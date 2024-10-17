#!/usr/bin/env bash

echo "sample_id,fastq_1,fastq_2" >fq.csv

for fqdir in $@; do
fqdir=$(realpath $fqdir)
	paste -d "," <(find $fqdir/ | egrep "_S[0-9]+(_L[0-9]+)?_R1_" | sort | sed "s/.*\///g; s/_S[0-9]\+_R1_.*//; s/_S[0-9]\+_L[0-9]\+_R1_.*//") <(find $fqdir/ | egrep "_S[0-9]+(_L[0-9]+)?_R1_" | sort ) <(find $fqdir/ | egrep "_S[0-9]+(_L[0-9]+)?_R2_" | sort )
done | egrep -vi "undetermined" >>fq.csv


