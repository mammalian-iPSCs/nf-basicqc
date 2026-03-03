/*
========================================================================================
    SORTMERNA Module
========================================================================================
    SortMeRNA - Tool for filtering ribosomal RNA from sequencing reads
    Adapted from nf-core/modules sortmerna module (v4.3.7)
    Used here for rRNA quantification only (log output for MultiQC, reads not saved)

    Note: SortMeRNA computes index file names from a hash of the reference file PATH,
    not its content. This means a pre-built index cannot be reused across different
    work directories (where staged FASTAs have different paths). Therefore each SORTMERNA
    job builds its own index in its work directory (--index 1). Nextflow's caching
    (-resume) avoids rebuilding on subsequent runs with the same inputs.
*/

process SORTMERNA {
    tag "$sample"
    label 'process_high'
    publishDir "${params.outdir}/sortmerna", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    tuple val(sample), path(reads)
    path(fastas)        // collected list of rRNA database FASTA files
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
        --index 1 \\
        ${paired_cmd} \\
        ${out2_cmd} \\
        $args

    mv rRNA_reads.log ${prefix}.sortmerna.log

    # Optionally save rRNA reads for downstream classification
    if [[ "${save_rrna}" == "true" ]]; then
        if [[ "${paired}" == "true" ]]; then
            mv rRNA_reads_fwd.fastq.gz ${prefix}_rrna_1.fastq.gz
            mv rRNA_reads_rev.fastq.gz ${prefix}_rrna_2.fastq.gz
        else
            mv rRNA_reads.fastq.gz ${prefix}_rrna.fastq.gz
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
