#!/bin/bash
#SBATCH --job-name=sortmerna_db_test
#SBATCH --output=sortmerna_db_test_%j.out
#SBATCH --error=sortmerna_db_test_%j.err
#SBATCH --partition=genD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --qos=normal
#SBATCH --time=4:00:00
#SBATCH --exclude=cnd75,cnd73

# Test each SortMeRNA rRNA database file individually on the CGLZOO_01 test samples.
# Useful for understanding which databases contribute rRNA signal across taxa.
# Usage: sbatch ribodetector_test/run_sortmerna_db_test.sh

DB_DIR="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/rRNA_indices"
FASTQ_DIR="/scratch_isilon/groups/compgen/data_transfer/CGLZOO_01/20241128/FASTQ"
OUTDIR="/scratch_isilon/groups/compgen/lwange/nf-basicqc/ribodetector_test/sortmerna_per_db"
SINGULARITY_CACHE="/scratch_isilon/groups/compgen/lwange/singularity/basicqc"
SORTMERNA_CONTAINER="https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/919d9c8f5f2c3221a94efe96b81bde0c953c13ebb0a1eca6b690b90666006cad/data"
SEQTK_CONTAINER="https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2"
SUBSAMPLE=100000

declare -A SAMPLE_LABELS=(
    ["HFYMJDSXC_1_8bp-UDP0032"]="BB1523 (Callithrix geoffroyi)"
    ["HFYMJDSXC_1_8bp-UDP0034"]="BB1525 (Gorilla gorilla)"
)
SAMPLES=("HFYMJDSXC_1_8bp-UDP0032" "HFYMJDSXC_1_8bp-UDP0034")

mkdir -p "${OUTDIR}/subsampled"

module load singularity 2>/dev/null || true

run_sortmerna() { singularity exec --bind /scratch_isilon "${SORTMERNA_CONTAINER}" "$@"; }
run_seqtk()    { singularity exec --bind /scratch_isilon "${SEQTK_CONTAINER}"    "$@"; }

# ── Step 1: subsample reads once ────────────────────────────────────────────
echo "$(date) Subsampling reads to ${SUBSAMPLE}..."
for sample in "${SAMPLES[@]}"; do
    sub_r1="${OUTDIR}/subsampled/${sample}_1.fastq.gz"
    sub_r2="${OUTDIR}/subsampled/${sample}_2.fastq.gz"
    if [ ! -f "${sub_r1}" ]; then
        seed=$(echo "${sample}" | cksum | cut -d' ' -f1)
        run_seqtk seqtk sample -s${seed} "${FASTQ_DIR}/${sample}_1.fastq.gz" ${SUBSAMPLE} \
            | gzip -c > "${sub_r1}"
        run_seqtk seqtk sample -s${seed} "${FASTQ_DIR}/${sample}_2.fastq.gz" ${SUBSAMPLE} \
            | gzip -c > "${sub_r2}"
        echo "  Subsampled ${sample}"
    else
        echo "  Skipping ${sample} (already subsampled)"
    fi
done

# ── Step 2: run SortMeRNA with each database individually ───────────────────
shopt -s nullglob
db_files=("${DB_DIR}"/*.fasta "${DB_DIR}"/*.fa "${DB_DIR}"/*.fna)
shopt -u nullglob

if [ ${#db_files[@]} -eq 0 ]; then
    echo "ERROR: no .fasta/.fa/.fna files found in ${DB_DIR}" >&2
    exit 1
fi

echo ""
echo "$(date) Found ${#db_files[@]} database file(s):"
for db in "${db_files[@]}"; do echo "  $(basename ${db})"; done
echo ""

for db_fasta in "${db_files[@]}"; do
    db_name=$(basename "${db_fasta}" | sed 's/\.[^.]*$//')
    db_outdir="${OUTDIR}/${db_name}"
    mkdir -p "${db_outdir}"

    echo "$(date) === ${db_name} ==="

    for sample in "${SAMPLES[@]}"; do
        workdir="${db_outdir}/${sample}_workdir"
        mkdir -p "${workdir}"
        log_file="${db_outdir}/${sample}.sortmerna.log"

        echo "  ${sample}..."
        run_sortmerna sortmerna \
            --ref "${db_fasta}" \
            --reads "${OUTDIR}/subsampled/${sample}_1.fastq.gz" \
            --reads "${OUTDIR}/subsampled/${sample}_2.fastq.gz" \
            --threads ${SLURM_CPUS_PER_TASK} \
            --workdir "${workdir}" \
            --fastx \
            --aligned "${workdir}/rRNA_reads" \
            --paired_in \
            --out2

        mv "${workdir}/rRNA_reads.log" "${log_file}"
        rm -rf "${workdir}"
    done
done

# ── Step 3: summary table ────────────────────────────────────────────────────
echo ""
echo "$(date) ============================================"
echo "         RESULTS: % rRNA per database"
echo "============================================"
printf "%-50s" "Database"
for sample in "${SAMPLES[@]}"; do
    printf "  %-30s" "${SAMPLE_LABELS[$sample]}"
done
echo ""
printf "%-50s" "--------"
for sample in "${SAMPLES[@]}"; do printf "  %-30s" "------------------------------"; done
echo ""

for db_fasta in "${db_files[@]}"; do
    db_name=$(basename "${db_fasta}" | sed 's/\.[^.]*$//')
    printf "%-50s" "${db_name}"
    for sample in "${SAMPLES[@]}"; do
        log="${OUTDIR}/${db_name}/${sample}.sortmerna.log"
        if [ -f "${log}" ]; then
            pct=$(grep -i "passing E-value threshold\|Total reads aligned" "${log}" \
                  | grep -oP '\([\d.]+%\)' | tr -d '()' | tail -1)
            printf "  %-30s" "${pct:-n/a}"
        else
            printf "  %-30s" "missing"
        fi
    done
    echo ""
done

echo ""
echo "$(date) Done. Logs in: ${OUTDIR}/"
