#!/bin/bash
#SBATCH --job-name=nf_basicqc_test
#SBATCH --output=nf_basicqc_test_%j.out
#SBATCH --error=nf_basicqc_test_%j.err
#SBATCH --partition=genD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --exclude=cnd75,cnd73

# BasicQC Nextflow Pipeline - Test Suite
# Submit with: sbatch submit_tests.sh
# For Kraken2 only: sbatch submit_tests.sh --kraken-only
# For fresh Kraken2 test (no resume): sbatch submit_tests.sh --kraken-fresh
# For FastQC only: sbatch submit_tests.sh --fastqc_only
# For full pipeline: sbatch submit_tests.sh --full
# For sex determination test: sbatch submit_tests.sh --sex
# For SortMeRNA rRNA test: sbatch submit_tests.sh --sortmerna

# Set paths
PIPELINE_DIR="/scratch_isilon/groups/compgen/lwange/nf-basicqc"
SLURM_CONFIG="/home/groups/compgen/lwange/isilon/lwange/singularity/basicqc/slurm.config"
FASTQ_SCREEN_CONF="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/FastQ_Screen_Genomes/FastQ_Screen_Genomes/fastq_screen.conf"
KRAKEN2_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_mtdna"
SEX_MARKERS_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/sex_markers/all_sex_markers.fasta"
SORTMERNA_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/rRNA_indices"

# Project metadata for MultiQC report
PROJECT_NAME="CGLZOO_01"
APPLICATION="RNA-seq"

cd $PIPELINE_DIR

echo "$(date) Starting BasicQC pipeline tests"
echo "========================================="

# Load nextflow if needed (uncomment/modify as needed for your cluster)
# module load nextflow
# module load singularity

arg1="${1:---full}"

# fastqc only pipeline test only
if [[ "$arg1" == "--fastqc_only" ]]; then
    echo "$(date) === Running full pipeline test ==="
    nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir test/results_fastqc \
        --skip_fastq_screen \
        --skip_kraken2 \
        --project_name $PROJECT_NAME \
        --application $APPLICATION \
        -profile singularity -c $SLURM_CONFIG \
        -resume
    echo "$(date) === FASTQC test complete ==="
    echo "Check results in: $PIPELINE_DIR/test/results_fastqc"
    exit 0
fi

# Check for --kraken-only flag
if [[ "$arg1" == "--kraken-only" ]]; then
    echo "$(date) === Running Kraken2 batch test only ==="
    nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir test/results_kraken2 \
        --skip_fastqc \
        --skip_fastq_screen \
        --kraken2_db $KRAKEN2_DB \
        --kraken2_subsample 100000 \
        --project_name $PROJECT_NAME \
        --application $APPLICATION \
        -profile singularity -c $SLURM_CONFIG \
        -resume
    echo "$(date) === Kraken2 test complete ==="
    echo "Check results in: $PIPELINE_DIR/test/results_kraken2"
    exit 0
fi

# Fresh Kraken2 test without resume (use after config changes)
if [[ "$arg1" == "--kraken-fresh" ]]; then
    echo "$(date) === Running fresh Kraken2 batch test (no resume) ==="
    # Clean previous work dirs to ensure fresh job submission
    echo "Cleaning previous work directories..."
    rm -rf work/
    nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir test/results_kraken2 \
        --skip_fastqc \
        --skip_fastq_screen \
        --kraken2_db $KRAKEN2_DB \
        --kraken2_subsample 100000 \
        --project_name $PROJECT_NAME \
        --application $APPLICATION \
        -profile singularity -c $SLURM_CONFIG
    echo "$(date) === Kraken2 test complete ==="
    echo "Check results in: $PIPELINE_DIR/test/results_kraken2"
    exit 0
fi

# Full pipeline test only
if [[ "$arg1" == "--full" ]]; then
    echo "$(date) === Running full pipeline test ==="
    nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir test/results_full \
        --fastq_screen_conf $FASTQ_SCREEN_CONF \
        --kraken2_db $KRAKEN2_DB \
        --kraken2_subsample 100000 \
        --sex_markers_db $SEX_MARKERS_DB \
        --sortmerna_db $SORTMERNA_DB \
        --sortmerna_subsample 100000 \
        --project_name $PROJECT_NAME \
        --application $APPLICATION \
        -profile singularity -c $SLURM_CONFIG \
        -resume
    echo "$(date) === Full pipeline test complete ==="
    echo "Check results in: $PIPELINE_DIR/test/results_full"
    exit 0
fi

# SortMeRNA rRNA quantification test only
if [[ "$arg1" == "--sortmerna" ]]; then
    echo "$(date) === Running SortMeRNA rRNA quantification test ==="
    nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir test/results_sortmerna \
        --skip_fastqc \
        --skip_fastq_screen \
        --skip_kraken2 \
        --skip_sex_determination \
        --sortmerna_db $SORTMERNA_DB \
        --sortmerna_subsample 100000 \
        --project_name $PROJECT_NAME \
        --application $APPLICATION \
        -profile singularity -c $SLURM_CONFIG \
        -resume
    echo "$(date) === SortMeRNA test complete ==="
    echo "Check results in: $PIPELINE_DIR/test/results_sortmerna"
    echo "Index saved to:   $PIPELINE_DIR/test/results_sortmerna/sortmerna/idx"
    echo "  (reuse with --sortmerna_index to skip rebuilding)"
    exit 0
fi

# Sex determination test only (with Kraken2)
if [[ "$arg1" == "--sex" ]]; then
    echo "$(date) === Running sex determination test ==="
    nextflow run main.nf \
        --input test/test_samplesheet.csv \
        --outdir test/results_sex \
        --skip_fastqc \
        --skip_fastq_screen \
        --kraken2_db $KRAKEN2_DB \
        --kraken2_subsample 100000 \
        --sex_markers_db $SEX_MARKERS_DB \
        --project_name $PROJECT_NAME \
        --application $APPLICATION \
        -profile singularity -c $SLURM_CONFIG \
        -resume
    echo "$(date) === Sex determination test complete ==="
    echo "Check results in: $PIPELINE_DIR/test/results_sex"
    exit 0
fi


