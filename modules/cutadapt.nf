process CUTADAPT {
    label "process_high"
    module = ['cutadapt/4.9-GCCcore-12.3.0']

    tag "CUTADAPT on ${meta.id}"

    publishDir "${params.outdir}/trimmed_fastq/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    
    output:
    tuple val(meta), val(out_prefix), path( "output/${read1}" ), path( "output/${read2}" ), emit: fq
    tuple val(meta), val(out_prefix), path( "${out_prefix}.cutadapt.json" ), emit: js

    script:
    def args = task.ext.args ?: ""
    def adapters = params.adapters ?: "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    """
    mkdir output
    cutadapt -j ${task.cpus} \
    $args $adapters \
    --json=${out_prefix}.cutadapt.json \
		--nextseq-trim=20 \
		-m 20 \
		--overlap 3 \
		-o output/$read1 \
    -p output/$read2 \
    $read1 $read2

    """
}
