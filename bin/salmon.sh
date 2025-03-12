#!/usr/bin/env bash

cpus=$1; shift
bam=$1; shift
tx_fa=$1; shift
out_prefix=$1; shift

salmon quant -t $tx_fa -l A -a $bam -p $cpus \
--numBootstraps 200 \
--gcBias --seqBias --reduceGCMemory --biasSpeedSamp 5 --incompatPrior 0.0 -o $out_prefix
