#!/usr/bin/env bash

sample_id=$1; shift
gtf=$1; shift
bam=$1; shift
read_type=$1; shift
strand=$1; shift
nthread=$1; shift

if [[ $read_type == PE ]]; then
	featureCounts -T $nthread -p -C --minOverlap 1 -a $gtf -F GTF -t exon -g gene_id -s $strand -o ${sample_id}.exonicFragmentCounts.txt $bam
else
    featureCounts -T $nthread --minOverlap 1 -a $gtf -F GTF -t exon -g gene_id -s $strand -o ${sample_id}.exonicReadCounts.txt $bam
fi

