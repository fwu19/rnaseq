params.outdir = 'Analysis/'

process MULTIQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['MultiQC/1.9-foss-2019b-Python-3.7.4']


    tag "MultiQC on all samples"

    publishDir "${outdir}/MultiQC/", mode: 'copy'

   input:
    tuple val(sample_id), path(bam, stageAs: 'STAR/*')
    path fastqc, stageAs: 'FastQC/*'
    path rseqc, stageAs: 'RSeQC/*'
    path rnaseqc, stageAs: 'RNA-SeQC/*'
    path outdir

    output:
    path '*'

    script:
    """
    multiqc ./ -o . --ignore _STARpass1/ -f 

    """
}