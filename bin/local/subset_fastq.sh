#!/usr/bin/env bash

outdir=$1; shift
for fq in "$@"; do
    zcat $fq | head -n 40000 | gzip -c >$outdir/$(basename $fq)
done