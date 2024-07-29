#!/usr/bin/env bash

sampleId=$1; shift
sampleBam=$1/${sampleId}.bam; shift
gtfQC=$1; shift
strand=${1:-2}; shift
readType=${1:-PE}; shift

if [[ $strand -eq 2 ]]; then
    moreArgs="--stranded=RF"
fi

if [[ $readType -eq "SE" ]]; then
	moreArgs="${moreArgs} -u"
fi

rnaseqc $gtfQC $sampleBam ./ -s $sampleId $moreArgs

	