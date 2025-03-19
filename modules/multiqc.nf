process MULTIQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['MultiQC/1.21-foss-2023a']


    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    path ( multiqc_custom_config)
    path ( 'star_log/*' )
    path ( 'star_count/*' )
    path ( 'fastqc/*' )
    path ( 'fastqc_trimmed/*' )
    path ( 'rseqc/*' )
    path ( 'rnaseqc/*' )
    path ( 'gatk/*' )

    output:
    path ('multiqc_data/*'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f --ignore _STARpass1/ $custom_config .
    """
}