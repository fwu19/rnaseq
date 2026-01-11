#!/usr/bin/env bash

bam=$1; shift
tx_bed=$1; shift
infer_experiment.py -i $bam -r $tx_bed $@ >out.txt

a=$(tail -n 2 out.txt| head -n 1| sed "s/.*: //g")
b=$(tail -n 1 out.txt| sed "s/.*: //g")

echo "strand,read_type" >infer_experiment.csv
if [[ $a > 0.6 ]]; then
    echo -n 1 >>infer_experiment.csv
elif [[ $b > 0.6 ]]; then
    echo -n 2 >>infer_experiment.csv
else
    echo -n 0 >>infer_experiment.csv
fi

type=$(tail -n 4 out.txt | head -n 1)
if [[ $type =~ 'PairEnd' ]]; then
    echo ',PE' >>infer_experiment.csv
else
    echo ',SE' >>infer_experiment.csv
fi
