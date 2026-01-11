
process MULTIQC {
    module = ['MultiQC/1.21-foss-2023a']

    label "process_medium"

    tag "MultiQC on all samples"

    input:
    path ( multiqc_config )
    path (fastqc) //( 'fastqc/*' )
    path (cutadapt) //( 'cutadapt/*' )
    path (fastp) //( 'fastp/*' )
    path (fastqc_trimmed) //( 'fastqc_trimmed/*' )
    path (rseqc) //( 'rseqc/*' )
    path (rnaseqc) //( 'rnaseqc/*' )
    path (gatk) //( 'gatk/*' )
    path (star_log) //( 'star_log/*' )
    path (star_log_host) //( 'star_log_host/*' )
    path (star_count) //( 'star_count/*' )
    path (samtools) //( 'samtools/*' )
    path (samtools_host) //( 'samtools_host/*' )
    path (samtools_xeno) //( 'samtools_xeno/*' )
    path (salmon) //( 'salmon/*' )

    output:
    path ('multiqc_data/', type: 'dir'), emit: data
    path ('multiqc_data/', type: 'dir' )
    path ('multiqc_report.html')
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f --ignore _STARpass1/ --config $multiqc_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    """
}