#!/usr/bin/env bash

: <<'END_COMMENT'
The input POSITIONAL parameters are as follows:
sample - sample_id
STARref - path/to/STAR/index
GTFfile - path/to/GTF/file
nthread - number of threads
reads - <path/to/read1.fastq.gz> [path/to/reads2.fastq.gz]
RG - define read group
END_COMMENT

sampleId=$1; shift
STARref=$1; shift
GTFfile=$1; shift
nthread=$1; shift

fq1=$(find fastq/ | egrep "merged_R1.fastq.gz|S[0-9]+(_L[0-9]+)?_R1_.*fastq.gz" )
fq2=$(find fastq/ | egrep "merged_R2.fastq.gz|S[0-9]+(_L[0-9]+)?_R2_.*fastq.gz" )
[[ -d $sampleId ]] || mkdir $sampleId
RG="ID:$sampleId LB:$sampleId SM:$sampleId PL:illumina PU:$sampleId CN:FredHutch"
STAR \
    --outFileNamePrefix "$sampleId/" \
	--genomeDir "${STARref}" \
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
	--sjdbGTFfile "${GTFfile}" \
	--outReadsUnmapped None # Could use "within", compatible with outSAMunmapped


samtools sort -T . -@ 6 -o $sampleId/Aligned.sortedByCoord.out.bam $sampleId/Aligned.out.bam
rm $sampleId/Aligned.out.bam

mv $sampleId/Aligned.sortedByCoord.out.bam $sampleId/${sampleId}.bam
samtools index $sampleId/${sampleId}.bam

