process WRITE_CSV {

    label 'process_single'

    tag "Write ${out_csv}."

    input:
    val (in_data)
    val (out_csv)

    output:
    path (out_csv)

    script:
    def header = in_data[0].keySet().join(',')
    def rows = in_data.collect { row ->
        row.values().join(',')
    }.join('\n')

    """
    echo "${header}" > ${out_csv}
    echo "${rows}" >> ${out_csv}
    """
}
