
process XENGSORT_CLASSIFY {
    label "process_high"
    container "docker://quay.io/biocontainers/xengsort:2.2.1--pyhdfd78af_0" 


    tag "run xengsort on ${out_prefix}"

    publishDir "${params.outdir}/fastq/${params.genome}_${params.genome_host}", mode: 'copy'
    
    input:
    tuple val(meta), val(out_prefix), path(read1), path(read2)
    path (xengsort_dir)

    output:
    tuple val(meta), val(out_prefix), path ("${out_prefix}-graft.1.fq.gz"), path ("${out_prefix}-graft.2.fq.gz"), emit: graft_fq
    tuple val(meta), val(out_prefix), path ("${out_prefix}-host.1.fq.gz"), path ("${out_prefix}-host.2.fq.gz"), emit: host_fq
    path( "${out_prefix}.xengsort.log" ), emit: log
    path("*")
    path ('versions.yml'), emit: versions
    
    script:
    def args = task.ext.args ?: ""
    """
    xengsort_index=\$(cd ${xengsort_dir}; ls *.hash | sed 's/[.]hash\$//')
    xengsort classify -T ${task.cpus} --index ${xengsort_dir}/\${xengsort_index} --fastq $read1 --pairs $read2 --prefix $out_prefix | tee ${out_prefix}.xengsort.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Xengsort: \$(xengsort --version 2>&1)
    END_VERSIONS
    
    """
}
