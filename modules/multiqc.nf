process MULTIQC {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['MultiQC/1.21-foss-2023a']


    tag "MultiQC on all samples"

    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    path ( multiqc_config)
    path ( 'fastqc/*' )
    path ( 'cutadapt/*' )
    path ( 'fastqc_trimmed/*' )
    path ( 'rseqc/*' )
    path ( 'rnaseqc/*' )
    path ( 'gatk/*' )
    path ( 'star_log/*' )
    path ( 'star_count/*' )
    path ( 'salmon/*' )

    output:
    path ('multiqc_data/*'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f --ignore _STARpass1/ --config $multiqc_config .
    """
}