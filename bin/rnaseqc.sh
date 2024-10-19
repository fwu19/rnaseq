#!/usr/bin/env bash

sampleId=$1; shift
gtfQC=$1; shift
strand=${1:-2}; shift
readType=${1:-PE}; shift

if [[ $strand -eq 2 ]]; then
    moreArgs="--stranded=RF"
fi

if [[ $readType -eq "SE" ]]; then
	moreArgs="${moreArgs} -u"
fi

for sampleBam in "$@"; do
    rnaseqc $gtfQC $sampleBam ./ -s $sampleId $moreArgs
done
	