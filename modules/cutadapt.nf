process CUTADAPT {
    label "process_high"
    module = ['cutadapt/4.9-GCCcore-12.3.0']

    tag "CUTADAPT on ${meta.id}"

    publishDir "${params.outdir}/trimmed_fastq/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    
    output:
    tuple val(meta), val(out_prefix), path("trimmed_fastq/${read1}"), path("trimmed_fastq/${read2}"), emit: fq
    tuple val(meta), val(out_prefix), path( "${out_prefix}.cutadapt.json" ), emit: js

    script:
    def args = task.ext.args ?: ""
    def adapter_list = params.adapters.split(',').collect()
    if ( adapter_list.size == 1){
        adapter1 = adapter_list[0]
        adapter2 = adapter_list[0]
    } else {
        adapter1 = adapter_list[0]
        adapter2 = adapter_list[1]
    }
    """
    mkdir trimmed_fastq
    cutadapt -j ${task.cpus} \
    $args \
    -a $adapter1 -A $adapter2 \
    --json=${out_prefix}.cutadapt.json \
		--nextseq-trim=20 \
		-m 20 \
		--overlap 3 \
		-o trimmed_fastq/$read1 \
    -p trimmed_fastq/$read2 \
    $read1 $read2

    """
}
