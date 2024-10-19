#!/usr/bin/env bash

module purge all
module load MultiQC/1.9-foss-2019b-Python-3.7.4
outdir=$1; shift
[[ -d $outdir ]] || mkdir -p $outdir
outdir=$( realpath $outdir )
multiqc $@ -o $outdir --ignore-symlinks --ignore _STARpass1/ -f