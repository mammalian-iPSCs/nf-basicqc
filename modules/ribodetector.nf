/*
========================================================================================
    RIBODETECTOR Module
========================================================================================
    RiboDetector - ML-based rRNA detection without a reference database
    Adapted from nf-core/modules ribodetector module (v0.3.2)
    Used here for rRNA quantification only (log output for MultiQC, reads not saved)
*/

process RIBODETECTOR {
    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/ribodetector", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(sample), path(reads)
    val(length)    // sequenced read length (bp)

    output:
    tuple val(sample), path("*.ribodetector.log"), emit: log
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: sample
    def paired = reads instanceof List && reads.size() == 2
    def input_cmd  = paired ? "-i ${reads[0]} ${reads[1]}" : "-i ${reads}"
    def output_cmd = paired \
        ? "-o ${prefix}.nonrna_1.fastq.gz ${prefix}.nonrna_2.fastq.gz" \
        : "-o ${prefix}.nonrna.fastq.gz"
    """
    ribodetector_cpu \\
        ${input_cmd} \\
        ${output_cmd} \\
        -l ${length} \\
        -t ${task.cpus} \\
        -e rrna \\
        --chunk_size 256 \\
        --log ${prefix}.ribodetector.log \\
        ${args}

    # Clean up read files — only the log is needed for QC
    rm -f ${prefix}.nonrna*.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: \$(ribodetector_cpu --version | sed 's/ribodetector //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: sample
    """
    touch ${prefix}.ribodetector.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector: 0.3.2
    END_VERSIONS
    """
}
