#!/usr/bin/env bash

: <<'END_COMMENT'
The input POSITIONAL parameters are as follows:
sample_id - sample ID
star_ref - path/to/STAR/index
gtf - path/to/GTF/file
nthread - number of threads
fq1 - <path/to/read1.fastq.gz> 
fq2 - [path/to/reads2.fastq.gz]
RG - define read group
END_COMMENT

sample_id=$1; shift
star_ref=$1; shift
gtf=$1; shift
nthread=$1; shift
fq1=$1; shift
fq2=$1; shift


RG="ID:$sample_id LB:$sample_id SM:$sample_id PL:illumina PU:$sample_id CN:FredHutch"
STAR \
	--genomeDir "${star_ref}" \
	--readFilesIn ${fq1} ${fq2} \
	--readFilesCommand zcat \
	--runThreadN ${nthread} \
	--outFilterMultimapScoreRange 1 \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNmax 10 \
	--alignIntronMax 500000 \
	--alignMatesGapMax 1000000 \
	--sjdbScore 2 \
	--alignSJDBoverhangMin 1 \
	--genomeLoad NoSharedMemory \
	--outFilterMatchNminOverLread 0.33 \
	--outFilterScoreMinOverLread 0.33 \
	--sjdbOverhang 100 \
	--twopassMode Basic \
	--outSAMstrandField intronMotif \
	--outSAMattributes NH HI NM MD AS XS \
	--outSAMunmapped Within \
	--outSAMtype BAM Unsorted \
	--outSAMheaderHD @HD VN:1.4 \
	--outSAMattrRGline "${RG}" \
	--quantMode TranscriptomeSAM GeneCounts \
	--sjdbGTFfile "${gtf}" \
	--outReadsUnmapped None # Could use "within", compatible with outSAMunmapped


samtools sort -T . -@ 6 -o ${sample_id}.bam Aligned.out.bam
rm Aligned.out.bam

samtools index ${sample_id}.bam

