#!/usr/bin/env nextflow

/*
========================================================================================
    BasicQC Pipeline
========================================================================================
    A Nextflow pipeline for basic QC of Illumina FASTQ files
    - FastQC: Quality control metrics
    - FastQ Screen: Species/contamination detection
    - Kraken2: Taxonomic classification
    - MultiQC: Aggregate reports
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Import modules
include { FASTQC                                          } from './modules/fastqc'
include { FASTQ_SCREEN                                    } from './modules/fastq_screen'
include { SEQTK_SUBSAMPLE                                 } from './modules/seqtk_subsample'
include { SEQTK_SUBSAMPLE as SEQTK_SUBSAMPLE_RRNA         } from './modules/seqtk_subsample'
include { KRAKEN2                                         } from './modules/kraken2'
include { SUMMARIZE_KRAKEN2                               } from './modules/summarize_kraken2'
include { SEX_DETERMINATION                               } from './modules/sex_determination'
include { SUMMARIZE_SEX                                   } from './modules/sex_determination'
include { SORTMERNA_INDEX                                 } from './modules/sortmerna'
include { SORTMERNA                                       } from './modules/sortmerna'
include { RIBODETECTOR                                    } from './modules/ribodetector'
include { MULTIQC                                         } from './modules/multiqc'
include { PREPARE_MULTIQC_CONFIG                          } from './modules/prepare_multiqc_config'
include { SUMMARIZE_RESULTS                               } from './modules/summarize_results'

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================
     BasicQC Pipeline v1.0
    =========================================
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Mandatory arguments:
      --input           Path to input samplesheet (CSV format)
      --outdir          Output directory for results

    Optional arguments:
      --fastq_screen_conf   Path to fastq_screen configuration file
      --kraken2_db          Path to Kraken2 database
      --kraken2_subsample   Number of reads to subsample for Kraken2 (default: 5000000)
      --sex_markers_db      Path to sex marker FASTA for sex determination
      --sortmerna_db        Path to directory containing rRNA FASTA database files
      --rrna_subsample      Number of reads to subsample for rRNA tools (default: 1000000)
      --read_length         Sequenced read length in bp for RiboDetector (default: 150)
      --skip_fastqc         Skip FastQC step
      --skip_fastq_screen   Skip FastQ Screen step
      --skip_kraken2        Skip Kraken2 step
      --skip_sex_determination  Skip sex determination step
      --skip_sortmerna      Skip SortMeRNA rRNA quantification step
      --project_name        Project name for MultiQC report header (e.g., 'CGLZOO_01')
      --application         Application type for MultiQC header (e.g., 'RNA-seq')
      -profile              Configuration profile (singularity, docker, conda)

    Samplesheet format (CSV):
      sample,fastq_1,fastq_2,sample_name,species
      HFYMJDSXC_1_8bp-UDP0032,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,BB1523,Callithrix geoffroyi
      HFYMJDSXC_1_8bp-UDP0034,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,BB1525,Gorilla gorilla

    Note: sample_name and species columns are optional but enable better MultiQC grouping
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check mandatory parameters
if (!params.input) {
    error "Please provide an input samplesheet with --input"
}

if (!params.outdir) {
    error "Please provide an output directory with --outdir"
}

// Check input file exists
input_file = file(params.input)
if (!input_file.exists()) {
    error "Input samplesheet not found: ${params.input}"
}

/*
========================================================================================
    INPUT CHANNEL
========================================================================================
*/

def parse_samplesheet(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2) : null

            if (!fastq_1.exists()) {
                error "FASTQ file not found: ${row.fastq_1}"
            }
            if (fastq_2 && !fastq_2.exists()) {
                error "FASTQ file not found: ${row.fastq_2}"
            }

            return fastq_2 ? tuple(sample, [fastq_1, fastq_2]) : tuple(sample, [fastq_1])
        }
}

// Parse samplesheet for metadata (sample_name, species)
def parse_samplesheet_metadata(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            [
                fli: row.sample,
                sample_name: row.sample_name ?: row.sample,
                species: row.species ?: ''
            ]
        }
        .collect()
}

// Parse samplesheet for Kraken2 - one FASTQ per sample_name only
// Groups by sample_name and takes only the first entry to reduce runtime
def parse_samplesheet_kraken2(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_name = row.sample_name ?: row.sample
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = row.fastq_2 ? file(row.fastq_2) : null
            def species = row.species ?: ''

            return fastq_2 ?
                tuple(sample_name, species, [fastq_1, fastq_2]) :
                tuple(sample_name, species, [fastq_1])
        }
        .groupTuple(by: [0, 1])  // Group by sample_name and species
        .map { sample_name, species, reads_list ->
            // Take only the first FASTQ pair for this sample
            tuple(sample_name, species, reads_list[0])
        }
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Parse samplesheet and create input channel
    ch_reads = parse_samplesheet(params.input)

    // Parse sample metadata for MultiQC config
    ch_sample_metadata = parse_samplesheet_metadata(params.input)

    // Initialize empty channels for MultiQC
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { it[1] })
    }

    //
    // MODULE: FastQ Screen
    //
    if (!params.skip_fastq_screen) {
        // Check for fastq_screen config
        if (!params.fastq_screen_conf) {
            log.warn "No FastQ Screen config provided (--fastq_screen_conf). Skipping FastQ Screen."
        } else {
            ch_fastq_screen_conf = file(params.fastq_screen_conf)
            FASTQ_SCREEN(ch_reads, ch_fastq_screen_conf)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_SCREEN.out.txt.map { it[1] })
        }
    }

    //
    // MODULE: Kraken2 (with subsampling) - ONE FASTQ PER SAMPLE
    // Only processes one FASTQ per sample_name to reduce runtime
    //
    if (!params.skip_kraken2) {
        if (!params.kraken2_db) {
            log.warn "No Kraken2 database provided (--kraken2_db). Skipping Kraken2."
        } else {
            ch_kraken2_db = file(params.kraken2_db)

            // Parse samplesheet for Kraken2 - one FASTQ per sample_name
            ch_kraken2_reads = parse_samplesheet_kraken2(params.input)

            // Prepare channel for subsampling: tuple(sample_name, reads)
            ch_kraken2_for_subsample = ch_kraken2_reads
                .map { sample_name, species, reads -> tuple(sample_name, reads) }

            // Subsample reads before Kraken2 for efficiency
            SEQTK_SUBSAMPLE(ch_kraken2_for_subsample, params.kraken2_subsample)

            // Run Kraken2 per sample
            KRAKEN2(SEQTK_SUBSAMPLE.out.reads, ch_kraken2_db)
            // Note: We don't pass raw Kraken reports to MultiQC to avoid the default plot
            // with unclassified reads. Instead, we use our custom SUMMARIZE_KRAKEN2 output.

            // Collect Kraken2 metadata for summarization (sample_name -> species mapping)
            ch_kraken2_metadata = ch_kraken2_reads
                .map { sample_name, species, reads ->
                    [sample_name: sample_name, species: species]
                }
                .collect()

            // Summarize Kraken2 results for MultiQC
            // Creates: general stats columns + modified reports without unclassified
            ch_kraken2_reports = KRAKEN2.out.report.map { it[1] }.collect()
            SUMMARIZE_KRAKEN2(ch_kraken2_reports, ch_kraken2_metadata)
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_KRAKEN2.out.summary)
            // Pass modified Kraken reports (without unclassified) to MultiQC for interactive plot
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_KRAKEN2.out.classified_reports.flatten())
        }
    }

    //
    // MODULE: Sex Determination
    // Uses the same subsampled reads as Kraken2
    //
    if (!params.skip_sex_determination) {
        if (!params.sex_markers_db) {
            log.warn "No sex markers database provided (--sex_markers_db). Skipping sex determination."
        } else {
            ch_sex_markers_db = file(params.sex_markers_db)

            // Parse samplesheet for sex determination - same as Kraken2
            ch_sex_reads = parse_samplesheet_kraken2(params.input)

            // Prepare channel for subsampling (if not already subsampled by Kraken2)
            ch_sex_for_subsample = ch_sex_reads
                .map { sample_name, species, reads -> tuple(sample_name, reads) }

            // Use same subsampled reads as Kraken2 if available, otherwise subsample
            if (!params.skip_kraken2 && params.kraken2_db) {
                // Reuse subsampled reads from Kraken2
                ch_sex_subsampled = SEQTK_SUBSAMPLE.out.reads
            } else {
                // Need to subsample separately
                SEQTK_SUBSAMPLE(ch_sex_for_subsample, params.kraken2_subsample)
                ch_sex_subsampled = SEQTK_SUBSAMPLE.out.reads
            }

            // Determine species class from samplesheet (mammal, bird, etc.)
            // For now, we use 'unknown' and let the module infer from marker hits
            ch_sex_with_class = ch_sex_subsampled
                .map { sample_name, reads -> tuple(sample_name, reads, 'unknown') }

            // Run sex determination
            SEX_DETERMINATION(
                ch_sex_subsampled,
                ch_sex_markers_db,
                'unknown'  // species class - inferred from markers
            )

            // Collect metadata for summarization
            ch_sex_metadata = ch_sex_reads
                .map { sample_name, species, reads ->
                    [sample_name: sample_name, species: species]
                }
                .collect()

            // Summarize sex determination results for MultiQC
            ch_sex_results = SEX_DETERMINATION.out.results.map { it[1] }.collect()
            SUMMARIZE_SEX(ch_sex_results, ch_sex_metadata)
            ch_multiqc_files = ch_multiqc_files.mix(SUMMARIZE_SEX.out.summary)
        }
    }

    //
    // rRNA quantification: SortMeRNA and/or RiboDetector
    // Subsampled reads are shared between both tools when both are enabled.
    //
    def run_sortmerna    = !params.skip_sortmerna && params.sortmerna_db
    def run_ribodetector = !params.skip_ribodetector

    if (run_sortmerna || run_ribodetector) {
        // Subsample once, shared by both tools
        SEQTK_SUBSAMPLE_RRNA(ch_reads, params.rrna_subsample)
        ch_rrna_reads = SEQTK_SUBSAMPLE_RRNA.out.reads
    }

    //
    // MODULE: SortMeRNA
    //
    if (run_sortmerna) {
        // Collect all rRNA FASTA files from the database directory
        ch_sortmerna_fastas = Channel
            .fromPath("${params.sortmerna_db}/*.{fasta,fa,fna}")
            .collect()

        // Build index once or reuse a pre-built one.
        // On first run the index is published to <outdir>/sortmerna/idx —
        // pass it as --sortmerna_index on subsequent runs to skip rebuilding.
        if (params.sortmerna_index) {
            ch_sortmerna_index = Channel.value(file(params.sortmerna_index))
        } else {
            SORTMERNA_INDEX(ch_sortmerna_fastas)
            ch_sortmerna_index = SORTMERNA_INDEX.out.index
        }

        SORTMERNA(ch_rrna_reads, ch_sortmerna_fastas, ch_sortmerna_index)
        ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log.map { it[1] })
    } else if (!params.skip_sortmerna) {
        log.warn "No SortMeRNA database provided (--sortmerna_db). Skipping SortMeRNA."
    }

    //
    // MODULE: RiboDetector
    //
    if (run_ribodetector) {
        RIBODETECTOR(ch_rrna_reads, params.read_length)
        ch_multiqc_files = ch_multiqc_files.mix(RIBODETECTOR.out.log.map { it[1] })
    }

    //
    // Generate MultiQC config with sample metadata
    //
    PREPARE_MULTIQC_CONFIG(
        ch_sample_metadata,
        params.project_name,
        params.application
    )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files
        .flatten()
        .collect()
        .filter { it.size() > 0 }
        .set { ch_multiqc_input }

    MULTIQC(
        ch_multiqc_input,
        PREPARE_MULTIQC_CONFIG.out.config
    )

    //
    // MODULE: Generate consolidated summary table
    //
    // Collect FastQC zip files
    ch_fastqc_for_summary = params.skip_fastqc
        ? Channel.of(file("NO_FASTQC"))
        : FASTQC.out.zip.map { it[1] }.collect()

    // Get Kraken2 summary (or placeholder)
    ch_kraken2_for_summary = (params.skip_kraken2 || !params.kraken2_db)
        ? Channel.of(file("NO_KRAKEN2"))
        : SUMMARIZE_KRAKEN2.out.summary

    // Get sex determination summary (or placeholder)
    ch_sex_for_summary = (params.skip_sex_determination || !params.sex_markers_db)
        ? Channel.of(file("NO_SEX"))
        : SUMMARIZE_SEX.out.summary

    // Parse sample info for summary
    ch_summary_sample_info = ch_sample_metadata.collect()

    SUMMARIZE_RESULTS(
        ch_fastqc_for_summary,
        ch_kraken2_for_summary,
        ch_sex_for_summary,
        ch_summary_sample_info
    )
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'failed'}"
    log.info "Results saved to: ${params.outdir}"
    log.info ""
    log.info "Key outputs:"
    log.info "  - Summary table: ${params.outdir}/summary/qc_summary.tsv"
    log.info "  - MultiQC report: ${params.outdir}/multiqc/"
    log.info ""
}
