process CUTADAPT {
    time = '1d'
    cpus = 6
    memory = '36G'
    module = ['cutadapt/4.9-GCCcore-12.3.0']


    tag "CUTADAPT on ${meta.id}"

    publishDir "${params.outdir}/trimmed_fastq/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2)
    

    output:
    tuple val(meta), path("output/${read1.name}"), path("output/${read2.name}"), emit: fq

    script:
    """
    mkdir output
    cutadapt -j ${task.cpus} \
		--nextseq-trim=20 \
			-m 20 \
				--overlap 3 \
					-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
						-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
							-o output/$read1 -p output/$read2 $read1 $read2

    """
}
