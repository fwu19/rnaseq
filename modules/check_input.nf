process CHECK_INPUT {
    input:
    path (file)
    
    output
    
    script:
    /*
    if (params.input == 'samples.csv')
        """
        echo "gzip -c $file > ${file}.gz"
        """
    else if (params.compress == 'bzip2')
        """
        echo "bzip2 -c $file > ${file}.bz2"
        """
    else
        throw new IllegalArgumentException("Unknown compressor $params.compress")
    */
    
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        // row is a list object
        .view { row -> "${row.sample_id}, ${row.fastq_1}, ${row.fastq_2}" }
        .set{reads_ch}

}

