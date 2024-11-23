#!/usr/bin/env bash

sample=$1; shift
bam=$1; shift
genomeFa=$1; shift
targetRegion=$1; shift

gatk --java-options -Xmx20g CollectHsMetrics \
      -I $bam \
      -O ${sample}.hs_metrics.txt \
      -R $genomeFa \
      -BAIT_INTERVALS $targetRegion \
      -TARGET_INTERVALS $targetRegion

