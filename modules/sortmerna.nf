/*
========================================================================================
    SORTMERNA Module
========================================================================================
    SortMeRNA - Tool for filtering ribosomal RNA from sequencing reads
    Adapted from nf-core/modules sortmerna module (v4.3.7)
    Used here for rRNA quantification (log output) and optionally to save rRNA reads
    for downstream Kraken2 classification (when save_rrna is true).

    Two processes:
      SORTMERNA_INDEX  — builds the index once per run from FASTA files; the output
                         is shared across all SORTMERNA jobs in the same run via
                         .first() (same pattern as nf-core/rnaseq).
      SORTMERNA        — runs per-sample classification using the shared index.

    Note: SortMeRNA hashes index files by the real (resolved) path of the reference
    FASTA. Within a single Nextflow run, all processes receive symlinks pointing to
    the same original files, so the hashes match and --index 0 works correctly.
    Cross-run index reuse is NOT supported — always let SORTMERNA_INDEX run fresh.
*/

process SORTMERNA_INDEX {
    tag "build_index"
    label 'process_high'

    input:
    path(fastas)

    output:
    path "idx"         , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fasta_list = fastas instanceof List ? fastas : [fastas]
    def refs_cmd   = fasta_list.collect { "--ref $it" }.join(' ')
    """
    sortmerna \\
        ${refs_cmd} \\
        --workdir . \\
        --index 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(sortmerna --version 2>&1 | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+" | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p idx
    touch idx/stub.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: 4.3.7
    END_VERSIONS
    """
}

process SORTMERNA {
    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/sortmerna", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(sample), path(reads)
    path(fastas)        // collected list of rRNA database FASTA files
    path(index)         // pre-built index directory (from SORTMERNA_INDEX)
    val(save_rrna)      // whether to save rRNA reads for downstream classification

    output:
    tuple val(sample), path("*.sortmerna.log"),              emit: log
    tuple val(sample), path("*_rrna*.fastq.gz"), optional: true, emit: rrna_reads
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: sample
    def paired     = reads instanceof List && reads.size() == 2
    def reads_cmd  = paired ? reads.collect { r -> "--reads $r" }.join(' ') : "--reads $reads"
    def paired_cmd = paired ? '--paired_in' : ''
    def out2_cmd   = paired ? '--out2' : ''
    def fasta_list = fastas instanceof List ? fastas : [fastas]
    def refs_cmd   = fasta_list.collect { "--ref $it" }.join(' ')
    """
    sortmerna \\
        ${refs_cmd} \\
        ${reads_cmd} \\
        --threads $task.cpus \\
        --workdir . \\
        --fastx \\
        --aligned rRNA_reads \\
        --idx-dir ${index} \\
        --index 0 \\
        ${paired_cmd} \\
        ${out2_cmd} \\
        $args

    mv rRNA_reads.log ${prefix}.sortmerna.log

    # Optionally save rRNA reads for downstream classification
    if [[ "${save_rrna}" == "true" ]]; then
        if [[ "${paired}" == "true" ]]; then
            mv rRNA_reads_fwd.f*.gz ${prefix}_rrna_1.fastq.gz
            mv rRNA_reads_rev.f*.gz ${prefix}_rrna_2.fastq.gz
        else
            mv rRNA_reads.f*.gz ${prefix}_rrna.fastq.gz
        fi
    else
        rm -f rRNA_reads*.fastq.gz rRNA_reads*.fq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(sortmerna --version 2>&1 | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+" | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: sample
    """
    touch ${prefix}.sortmerna.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: 4.3.7
    END_VERSIONS
    """
}
