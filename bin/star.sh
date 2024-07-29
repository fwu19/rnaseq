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

sample=$1; shift
STARref=$(realpath $1); shift
GTFfile=$(realpath $1); shift
nthread=$1; shift
reads="$@"

[[ -d $sample ]] || mkdir $sample
RG="ID:$sample LB:$sample SM:$sample PL:illumina PU:$sample CN:FredHutch"
STAR \
    --outFileNamePrefix "$sample/" \
	--genomeDir "${STARref}" \
	--readFilesIn ${reads} \
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

 
samtools sort -T . -@ 6 -o $sample/Aligned.sortedByCoord.out.bam $sample/Aligned.out.bam
rm $sample/Aligned.out.bam
 
mv $sample/Aligned.sortedByCoord.out.bam $sample/${sample}.bam
samtools index $sample/${sample}.bam

