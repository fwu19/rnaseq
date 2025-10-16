#!/usr/bin/env bash

bam=$1; shift
tx_bed=$1; shift
infer_experiment.py -i $bam -r $tx_bed $@ >infer.txt

a=$(tail -n 2 infer.txt| head -n 1| sed "s/.*: //g")
b=$(tail -n 1 infer.txt| sed "s/.*: //g")

if [[ $a > 0.6 ]]; then
    echo 1 > strand.txt
elif [[ $b > 0.6 ]]; then
    echo 2 > strand.txt
else
    echo 0 > strand.txt
fi

type=$(tail -n 4 infer.txt | head -n 1)
if [[ $type =~ 'PairEnd' ]]; then
    echo 'PE' > read_type.txt
else
    echo 'SE' > read_type.txt
fi
