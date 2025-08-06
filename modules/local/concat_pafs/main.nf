process CONCAT_PAFS {
    tag "concat pafs"

    input:
    path pafs

    output:
    path "annotated.paf.tsv", emit: cat

    script:
    """
    echo "Received for concat:" >&2
    ls -lh ${pafs} >&2 || echo "No files passed!" >&2

    cat ${pafs} > annotated.paf.tsv
    """
}
