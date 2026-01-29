process WRITE_JSON {

    label 'process_single'

    tag "Write ${out_json}."

    input:
    val (in_data)
    val (out_json)

    output:
    path (out_json)

    script:
    def json = groovy.json.JsonOutput.prettyPrint(
                groovy.json.JsonOutput.toJson(in_data)
            )
    """
    cat << 'EOF' > ${out_json}
    ${json}
    EOF
    """
}
