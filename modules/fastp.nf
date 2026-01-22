process FASTP {
    label "process_high"

    container = 'quay.io/biocontainers/fastp:0.23.4--h125f33a_4'
    //module = ['fastp/0.23.4-GCC-13.2.0']

    tag "FASTP on ${meta.id}"

    publishDir "${params.outdir}/trimmed_fastq/", pattern: '*.gz', mode: 'copy'

    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    

    output:
    tuple val(meta), val(out_prefix), path("trimmed_fastq/${read1}"), path("trimmed_fastq/${read2}"), emit: fq
    tuple val(meta), val(out_prefix), path( "${out_prefix}.fastp.json" ), emit: js
    tuple val(meta), val(out_prefix), path( "${out_prefix}.fastp.html" ), emit: html
    path ('versions.yml'), emit: versions

    script:
    def args = task.ext.args ?: ""
    if (params.adapters){
        adapter_list = params.adapters.split(',').collect()
        if ( adapter_list.size == 1){
            adapter1 = adapter_list[0]
            adapter2 = adapter_list[0]
        } else {
            adapter1 = adapter_list[0]
            adapter2 = adapter_list[1]
        }

        """
        mkdir trimmed_fastq
        fastp -w ${task.cpus} \
        $args \
        --adapter_sequence $adapter1 --adapter_sequence_r2 $adapter2 \
        --in1 $read1 --in2 $read2 \
        --out1 trimmed_fastq/$read1 --out2 trimmed_fastq/$read2 \
        -j ${out_prefix}.fastp.json -h ${out_prefix}.fastp.html \
        --detect_adapter_for_pe -l 20 -g

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$( fastp --version 2>&1 | head -n 1 | sed -e "s/.* //g" )
        END_VERSIONS

        """

    }else{

        """
        mkdir trimmed_fastq
        fastp -w ${task.cpus} \
        $args \
        --in1 $read1 --in2 $read2 \
        --out1 trimmed_fastq/$read1 --out2 trimmed_fastq/$read2 \
        -j ${out_prefix}.fastp.json -h ${out_prefix}.fastp.html \
        --detect_adapter_for_pe -l 20 -g

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$( fastp --version 2>&1 | head -n 1 | sed -e "s/.* //g" )
        END_VERSIONS

        """
    }
}
