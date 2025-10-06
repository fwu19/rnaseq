#!/usr/bin/env bash

bam=$1; shift
tx_bed=$1; shift
infer_experiment.py -i $bam -r $tx_bed $@ >infer.txt
a=$(tail -n 2 infer.txt| head -n 1| sed "s/.*: //g")
b=$(tail -n 1 infer.txt| sed "s/.*: //g")

if [[ $a > 0.6 ]]; then
    echo 1
elif [[ $b > 0.6 ]]; then
    echo 2
else
    echo 0
fi