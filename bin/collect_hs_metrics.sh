#!/usr/bin/env bash

sample=$1; shift
bam=$1; shift
targetRegion=$1; shift

samtools faidx genome.fa
gatk --java-options -Xmx20g CreateSequenceDictionary -R genome.fa -O genome.dict

(samtools view -H $bam | sed "s/^@RG.*/@RG\tID:${sample}\tSM:${sample}/"; samtools view $bam | sed "s/RG:Z.*/RG:Z:${sample}/") | \
gatk --java-options -Xmx20g CollectHsMetrics \
      -I /dev/stdin/ \
      -O ${sample}.hs_metrics.txt \
      -R genome.fa \
      -BAIT_INTERVALS $targetRegion \
      -TARGET_INTERVALS $targetRegion

