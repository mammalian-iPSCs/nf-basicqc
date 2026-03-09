#!/bin/bash
#SBATCH --job-name=nf_basicqc
#SBATCH --output=nf_basicqc_%j.out
#SBATCH --error=nf_basicqc_%j.err
#SBATCH --partition=genD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --qos=marathon
#SBATCH --time=48:00:00
#SBATCH --exclude=cnd75,cnd73

# BasicQC Nextflow Pipeline - Production Run
#
# Usage:
#   sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name] [application]
#
# Example:
#   sbatch submit_pipeline.sh inputs/CGLZOO_01.csv results/CGLZOO_01 CGLZOO_01 "RNA-seq"
#
# Features:
#   - FastQC: Read quality metrics
#   - FastQ Screen: Contamination screening
#   - Kraken2: mtDNA-based species identification
#   - Sex determination: Genetic sex inference from sex-specific markers
#   - MultiQC: Aggregated report with General Stats table

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name] [application]"
    echo ""
    echo "Arguments:"
    echo "  samplesheet.csv  - CSV file with columns: sample,fastq_1,fastq_2,sample_name,species"
    echo "  output_dir       - Directory for pipeline outputs"
    echo "  project_name     - (optional) Project name for MultiQC report title"
    echo "  application      - (optional) Application type (e.g., 'RNA-seq', 'WGS')"
    exit 1
fi

SAMPLESHEET="$(realpath $1)"
OUTDIR="$(realpath -m $2)"
PROJECT_NAME="${3:-BasicQC}"
APPLICATION="${4:-Quality Control}"

# Set paths
PIPELINE_DIR="/scratch_isilon/groups/compgen/lwange/nf-basicqc"
SLURM_CONFIG="/home/groups/compgen/lwange/isilon/lwange/singularity/basicqc/slurm.config"
FASTQ_SCREEN_CONF="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/FastQ_Screen_Genomes/FastQ_Screen_Genomes/fastq_screen.conf"
KRAKEN2_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_mtdna"
SEX_MARKERS_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/sex_markers/all_sex_markers.fasta"
SORTMERNA_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/rRNA_indices/smr_v4.3_default_db.fasta"
RRNA_KRAKEN2_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_animal_rrna"

# Create run directory and work from there to avoid lock conflicts
RUN_DIR="$(dirname $OUTDIR)"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

# Build a valid Nextflow run name: project_application_MMDDHHMI (unique per launch)
RUN_BASE=$(echo "${PROJECT_NAME}_${APPLICATION}" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9]/_/g; s/__*/_/g; s/^_//; s/_$//')
RUN_NAME="${RUN_BASE}_$(date +%m%d%H%M)"

echo "$(date) Starting BasicQC pipeline"
echo "=================================="
echo "Samplesheet:  $SAMPLESHEET"
echo "Output dir:   $OUTDIR"
echo "Run dir:      $RUN_DIR"
echo "Project:      $PROJECT_NAME"
echo "Application:  $APPLICATION"
echo ""
echo "Databases:"
echo "  FastQ Screen:  $FASTQ_SCREEN_CONF"
echo "  Kraken2 mtDNA: $KRAKEN2_DB"
echo "  Sex markers:   $SEX_MARKERS_DB"
echo "  SortMeRNA:     $SORTMERNA_DB"
echo "  rRNA Kraken2:  $RRNA_KRAKEN2_DB"
echo ""

# Load modules if needed (uncomment/modify as needed)
# module load nextflow
# module load singularity

# Run the pipeline
nextflow run ${PIPELINE_DIR}/main.nf \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --fastq_screen_conf "$FASTQ_SCREEN_CONF" \
    --kraken2_db "$KRAKEN2_DB" \
    --kraken2_subsample 1000000 \
    --sex_markers_db "$SEX_MARKERS_DB" \
    --sortmerna_db "$SORTMERNA_DB" \
    --rrna_kraken2_db "$RRNA_KRAKEN2_DB" \
    --project_name "$PROJECT_NAME" \
    --application "$APPLICATION" \
    -name "$RUN_NAME" \
    -w "${RUN_DIR}/work" \
    -profile singularity \
    -c $SLURM_CONFIG \
    -resume

echo ""
echo "$(date) Pipeline complete"
echo "Results in: $OUTDIR"
echo ""
echo "Key outputs:"
echo "  - Summary table:  $OUTDIR/summary/qc_summary.tsv"
echo "  - MultiQC report: $OUTDIR/multiqc/${PROJECT_NAME}_multiqc_report.html"
echo "  - FastQC results: $OUTDIR/fastqc/"
echo "  - Kraken2 reports: $OUTDIR/kraken2/"
