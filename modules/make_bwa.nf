
process MAKE_BWA {
    label "process_high"
    module = ['BWA/0.7.17-GCCcore-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "Make BWA index"

    publishDir "${params.outdir}/references/$genome/BWA/", mode: 'copy'

    input:
    path ("genome/*") 
    
    output:
    path( "BWAIndex/", emit: dir)
    path( "BWAIndex/", type: 'dir')
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    cat genome/* > genome.fa
    mkdir BWAIndex
    bwa index $args -p BWAIndex/genome genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BWA: \$(BWA --version | head -n 1)
    END_VERSIONS

    """
}
