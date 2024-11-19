process MULTIQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['MultiQC/1.9-foss-2019b-Python-3.7.4']


    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    tuple path ("*", stageAs: 'star/*' )
    tuple path ("*", stageAs: 'fastqc/*' )
    tuple path ("*", stageAs: 'rseqc/*' )
    tuple path ("*", stageAs: 'rnaseqc/*' )
    tuple path ("*", stageAs: 'gatk/*' )

    output:
    path ('multiqc_data/*')
    path ('multiqc_report.html')

    script:
    """
    multiqc -o . --ignore _STARpass1/ -f */

    """
}