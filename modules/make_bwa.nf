
process MAKE_BWA {
    time = '1d'
    cpus = 12
    memory = '48G'
    module = ['BWA/0.7.17-GCCcore-12.2.0', 'SAMtools/1.17-GCC-12.2.0']


    tag "Make BWA index"

    publishDir "${params.outdir}/references/$genome/BWA/", mode: 'copy'

    input:
    path ("genome/*") 
    
    output:
    path( "BWA/", emit: dir)
    path( "BWA/", type: 'dir')
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    """
    cat genome/* > genome.fa
    mkdir BWA
    bwa index $args -p BWA/genome genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        BWA: \$(BWA --version | head -n 1)
    END_VERSIONS

    """
}
