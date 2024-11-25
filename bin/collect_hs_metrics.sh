#!/usr/bin/env bash

sample=$1; shift
bam=$1; shift
genomeFa=$1; shift
targetRegion=$1; shift

(samtools view -H $bam | sed "s/^@RG.*/@RG\tID:${sample}\tSM:${sample}/"; samtools view $bam | sed "s/RG:Z.*/RG:Z:${sample}/") >${sample}.bam
samtools index ${sample}.bam

gatk --java-options -Xmx20g CollectHsMetrics \
      -I ${sample}.bam \
      -O ${sample}.hs_metrics.txt \
      -R $genomeFa \
      -BAIT_INTERVALS $targetRegion \
      -TARGET_INTERVALS $targetRegion

