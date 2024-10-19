process MULTIQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['MultiQC/1.9-foss-2019b-Python-3.7.4']


    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    tuple path ("*", stageAs: 'STAR/*')
    tuple path ("*", stageAs: 'FastQC/*')
    tuple path ("*", stageAs: 'RSeQC/*')
    tuple path ("*", stageAs: 'RNA-SeQC/*')

    output:
    path ('multiqc_data/*')
    path ('multiqc_report.html')

    script:
    """
    multiqc */ -o . --ignore _STARpass1/ -f 

    """
}