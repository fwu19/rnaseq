params.reads = "fastq/*_S[0-9]*_{R1,R2}*fastq.gz"
star = "/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/STAR2Index"
gtf = "/fh/fast/_SR/Genomics/proj/fwu/reference/hg38/genes/gencode.v38.annotation.proteinCoding_lncRNA.gtf"
params.outdir = 'Analysis'

process STAR {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['STAR/2.7.7a-GCC-10.2.0', 'SAMtools/1.11-GCC-10.2.0']


    tag "STAR on $sample_id"

    publishDir "${outdir}/STAR/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path(star)
    path(gtf)
    path(outdir)
    
    output:
    tuple val(sample_id), path("${sample_id}")

    script:
    """
    star.sh ${sample_id} ${star} ${gtf} $task.cpus ${reads1} ${reads2} 
    """
}
