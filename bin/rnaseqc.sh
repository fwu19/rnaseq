#!/usr/bin/env bash

sampleId=$1; shift
gtfQC=$1; shift
strand=$(tail -n 1 $1 | cut -f 1 -d , )
readType=$(tail -n 1 $1 | cut -f 2 -d , )

if [[ $strand -eq 2 ]]; then
    moreArgs="--stranded=RF"
elif [[ $strand -eq 2 ]]; then
    moreArgs="--stranded=FR"
fi

if [[ $readType -eq "SE" ]]; then
	moreArgs="${moreArgs} -u"
fi

for sampleBam in "$@"; do
    rnaseqc $gtfQC $sampleBam ./ -s $sampleId $moreArgs
done
	