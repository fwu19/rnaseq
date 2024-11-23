process MULTIQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['MultiQC/1.9-foss-2019b-Python-3.7.4']


    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    tuple path ( 'star/*' )
    tuple path ( 'fastqc/*' )
    tuple path ( 'rseqc/*' )
    tuple path ( 'rnaseqc/*' )
    tuple path ( 'gatk/*' )

    output:
    path ('multiqc_data/*'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f --ignore _STARpass1/ .
    """
}