/*
========================================================================================
    SEQTK_SUBSAMPLE Module
========================================================================================
    Subsample reads using seqtk for faster downstream analysis
*/

process SEQTK_SUBSAMPLE {
    tag "$sample"
    label 'process_low'

    input:
    tuple val(sample), path(reads)
    val(n_reads)

    output:
    tuple val(sample), path("*.subsampled.fastq.gz"), emit: reads
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: sample
    // Set seed for reproducibility (use sample hash for consistency across runs)
    def seed = sample.hashCode().abs() % 10000

    if (reads instanceof List && reads.size() == 2) {
        // Paired-end
        """
        seqtk sample -s${seed} ${reads[0]} ${n_reads} | gzip -c > ${prefix}_1.subsampled.fastq.gz \\
            || cp ${reads[0]} ${prefix}_1.subsampled.fastq.gz
        seqtk sample -s${seed} ${reads[1]} ${n_reads} | gzip -c > ${prefix}_2.subsampled.fastq.gz \\
            || cp ${reads[1]} ${prefix}_2.subsampled.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(seqtk 2>&1 | grep Version | sed 's/Version: //')
        END_VERSIONS
        """
    } else {
        // Single-end
        """
        seqtk sample -s${seed} ${reads} ${n_reads} | gzip -c > ${prefix}.subsampled.fastq.gz \\
            || cp ${reads} ${prefix}.subsampled.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(seqtk 2>&1 | grep Version | sed 's/Version: //')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: sample
    """
    touch ${prefix}_1.subsampled.fastq.gz
    touch ${prefix}_2.subsampled.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: 1.4
    END_VERSIONS
    """
}
