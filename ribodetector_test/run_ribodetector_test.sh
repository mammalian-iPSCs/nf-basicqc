#!/bin/bash
#SBATCH --job-name=ribodetector_test
#SBATCH --output=ribodetector_test_%j.out
#SBATCH --error=ribodetector_test_%j.err
#SBATCH --partition=genD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --threads-per-core=1
#SBATCH --mem=32G
#SBATCH --qos=normal
#SBATCH --time=2:00:00
#SBATCH --exclude=cnd75,cnd73

# RiboDetector CPU test on CGLZOO_01 test samples
# Usage: sbatch test/run_ribodetector_test.sh

OUTDIR="/scratch_isilon/groups/compgen/lwange/nf-basicqc/ribodetector_test/"
FASTQ_DIR="/scratch_isilon/groups/compgen/data_transfer/CGLZOO_01/20241128/FASTQ"

mkdir -p "${OUTDIR}"

echo "$(date) Starting RiboDetector test"
echo "Output: ${OUTDIR}"
echo ""


eval "$(conda shell.bash hook)"
conda activate ribodetector

# Sample 1: HFYMJDSXC_1_8bp-UDP0032 (BB1523 - Callithrix geoffroyi)
echo "$(date) Processing UDP0032 (Callithrix geoffroyi)..."
ribodetector_cpu \
    -t 20 \
    -l 150 \
    -i "${FASTQ_DIR}/HFYMJDSXC_1_8bp-UDP0032_1.fastq.gz" \
       "${FASTQ_DIR}/HFYMJDSXC_1_8bp-UDP0032_2.fastq.gz" \
    -e rrna \
    --chunk_size 256 \
    -o "${OUTDIR}/UDP0032.nonrrna_1.fq" \
       "${OUTDIR}/UDP0032.nonrrna_2.fq" \
    --rrna "${OUTDIR}/UDP0032.rrna_1.fq" \
            "${OUTDIR}/UDP0032.rrna_2.fq"
echo "$(date) UDP0032 done"

# Sample 2: HFYMJDSXC_1_8bp-UDP0034 (BB1525 - Gorilla gorilla)
echo "$(date) Processing UDP0034 (Gorilla gorilla)..."
ribodetector_cpu \
    -t 20 \
    -l 150 \
    -i "${FASTQ_DIR}/HFYMJDSXC_1_8bp-UDP0034_1.fastq.gz" \
       "${FASTQ_DIR}/HFYMJDSXC_1_8bp-UDP0034_2.fastq.gz" \
    -e rrna \
    --chunk_size 256 \
    -o "${OUTDIR}/UDP0034.nonrrna_1.fq" \
       "${OUTDIR}/UDP0034.nonrrna_2.fq" \
    --rrna "${OUTDIR}/UDP0034.rrna_1.fq" \
            "${OUTDIR}/UDP0034.rrna_2.fq"
echo "$(date) UDP0034 done"

echo ""
echo "$(date) RiboDetector test complete"
echo ""
echo "To check % rRNA per sample, run:"
echo "  for s in UDP0032 UDP0034; do"
echo "    total=\$(cat \${OUTDIR}/\${s}.rrna_1.fq \${OUTDIR}/\${s}.nonrrna_1.fq | awk 'NR%4==1' | wc -l)"
echo "    rrna=\$(awk 'NR%4==1' \${OUTDIR}/\${s}.rrna_1.fq | wc -l)"
echo "    echo \"\$s: \$rrna / \$total rRNA reads\""
echo "  done"
