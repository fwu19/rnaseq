#!/usr/bin/env bash

cpus=$1; shift
genome_dir=$1; shift
gene_dir=$1; shift
overhang=$1; shift

cat $gene_dir/* 2>/dev/null | egrep "^#" > genes.gtf
cat $gene_dir/* 2>/dev/null | egrep -v "^#" >> genes.gtf
[[ -s genes.gtf ]] && moreArgs="--sjdbGTFfile genes.gtf --sjdbOverhang $overhang"

STAR --runThreadN ${cpus} \
	--runMode genomeGenerate \
	--genomeDir STAR2Index/ \
	--genomeFastaFiles $genome_dir/* \
	$moreArgs \
    $@
