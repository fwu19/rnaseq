process MAKE_SYMLINKS {

    label 'process_single'

    tag "make symlinks of params.srcdir to params.outdir."

    input:
    path(rel_f)   // rel_f is relative path, e.g. sub1/sub2/file.bam
    val(srcdir)
    val(outdir)

    """
    abs_f=\$(realpath $rel_f)
    mkdir -p "${outdir}/\$(dirname "$abs_f")"
    # make symlink that points back to the original file under src_dir
    ln -sf "$abs_f" "${outdir}/$rel_f"
    """
}
