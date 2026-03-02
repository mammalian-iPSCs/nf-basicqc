/*
========================================================================================
    SORTMERNA Module
========================================================================================
    SortMeRNA - Tool for filtering ribosomal RNA from sequencing reads
    Adapted from nf-core/modules sortmerna module (v4.3.7)
    Used here for rRNA quantification only (log output for MultiQC, reads not saved)

    Two processes:
      SORTMERNA_INDEX  — builds the index once from FASTA files; index is saved to
                         outdir and can be reused across runs via --sortmerna_index
      SORTMERNA        — runs per-sample classification using the pre-built index
*/

process SORTMERNA_INDEX {
    tag "build_index"
    label 'process_high'
    publishDir "${params.outdir}/sortmerna/idx", mode: 'copy'

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
    path(fastas)   // collected list of rRNA database FASTA files
    path(index)    // pre-built index directory

    output:
    tuple val(sample), path("*.sortmerna.log"), emit: log
    path "versions.yml"                        , emit: versions

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

    # Clean up read files — only the log is needed for QC
    rm -f rRNA_reads*.fastq.gz rRNA_reads*.fq.gz

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
