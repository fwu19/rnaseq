
process MAKE_XENGSORT {
    label "process_medium"
    container "docker://quay.io/biocontainers/xengsort:2.2.1--pyhdfd78af_0" 


    tag "Build Xengsort index for ${params.genome}_${params.genome_host}"

    publishDir "${params.outdir}/references/${params.genome}_${params.genome_host}/XengsortIndex", mode: 'copy'
    
    input:
    path("graft_genome/*")
    path("host_genome/*")

    output:
    path("XengsortIndex"), emit: dir
    path ('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-n 4_500_000_000 -k 25"
    """
    [[ -d XengsortIndex ]] && rm -rf XengsortIndex
    mkdir -p XengsortIndex
    cat graft_genome/* > graft.fa
    cat host_genome/* > host.fa
    xengsort index --index XengsortIndex/graft_host_genomes -H host.fa -G graft.fa $args | tee XengsortIndex/xengsort_index.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Xengsort: \$(xengsort --version 2>&1)
    END_VERSIONS

    
    """
}
