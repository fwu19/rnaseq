process FASTP {
    label "process_high"
    module = ['fastp/0.23.4-GCC-13.2.0']


    tag "FASTP on ${meta.id}"

    publishDir "${params.outdir}/trimmed_fastq/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    

    output:
    tuple val(meta), val(out_prefix), path("output/${read1}"), path("output/${read2}"), emit: fq
    tuple val(meta), val(out_prefix), path( "${out_prefix}.fastp.json" ), emit: js
    tuple val(meta), val(out_prefix), path( "${out_prefix}.fastp.html" ), emit: html

    script:
    def args = task.ext.args ?: ""
    def adapters = params.adapters ?: ""
    """
    mkdir output
    fastp -w ${task.cpus} \
    $args $adapters \
    --in1 $read1 --in2 $read2 \
    --out1 output/$read1 --out2 output/$read2 \
    -j ${out_prefix}.fastp.json -h ${out_prefix}.fastp.html \
    --detect_adapter_for_pe -l 20 -g

    """
}
