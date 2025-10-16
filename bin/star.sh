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
ext_rg="$1"; shift


[[ -d $sample_id ]] || mkdir -p $sample_id
cd $sample_id 

RG="ID:$sample_id SM:$sample_id LB:$sample_id PL:illumina PU:$sample_id $ext_rg"
if [[ -s ../$gtf ]]; then
	more_args="--quantMode TranscriptomeSAM GeneCounts 	--sjdbGTFfile ../${gtf}"
else
	more_args=""
fi
STAR $@ $more_args \
	--genomeLoad NoSharedMemory \
	--genomeDir "../${star_ref}" \
	--readFilesIn ../${fq1} ../${fq2} \
	--readFilesCommand zcat \
	--limitBAMsortRAM 1250000000 \
	--outReadsUnmapped None \
	--outSAMtype BAM Unsorted \
	--outSAMstrandField intronMotif \
	--outSAMunmapped Within \
	--outSAMattrRGline "${RG}" \
	--outFilterMultimapNmax 20 \
	--outFilterMultimapScoreRange 1 \
	--outFilterMatchNminOverLread 0.33 \
	--outFilterScoreMinOverLread 0.33 \
	--outFilterMismatchNmax 10 \
	--alignIntronMax 500000 \
	--alignMatesGapMax 1000000 \
	--alignSJDBoverhangMin 1 \
	--sjdbScore 2 \
	--twopassMode Basic \
	--runThreadN ${nthread} 

samtools sort -T . -@ ${nthread} -o Aligned.sortedByCoord.out.bam Aligned.out.bam
[[ -f Aligned.sortedByCoord.out.bam ]] && rm Aligned.out.bam
[[ -f Aligned.sortedByCoord.out.bam ]] && samtools index Aligned.sortedByCoord.out.bam

