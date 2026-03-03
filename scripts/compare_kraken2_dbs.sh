#!/bin/bash
#SBATCH --job-name=kraken2_compare
#SBATCH --cpus-per-task=8
#SBATCH --mem=350G
#SBATCH --time=8:00:00
#SBATCH --output=kraken2_compare_%j.log

# Compare Kraken2 results between full database and mtDNA-only database
# Runs both databases on the same subsampled reads and compares species assignments

set -e

# Configuration
FULL_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/kraken"
MTDNA_DB="/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/kraken/k2_mtdna"
OUTDIR="${1:-kraken2_comparison}"
SAMPLESHEET="${2:-test_samplesheet.csv}"
SUBSAMPLE_SIZE="${3:-100000}"

# Singularity containers
KRAKEN_CONTAINER="https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data"
SEQTK_CONTAINER="https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2"
SINGULARITY_CACHE="/scratch_isilon/groups/compgen/lwange/singularity/basicqc"

module load singularity 2>/dev/null || true
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"

run_kraken2() {
    singularity exec --bind /scratch_isilon "${KRAKEN_CONTAINER}" "$@"
}

run_seqtk() {
    singularity exec --bind /scratch_isilon "${SEQTK_CONTAINER}" "$@"
}

echo "=== Kraken2 Database Comparison ==="
echo "Full DB: ${FULL_DB}"
echo "mtDNA DB: ${MTDNA_DB}"
echo "Samplesheet: ${SAMPLESHEET}"
echo "Subsample size: ${SUBSAMPLE_SIZE}"
echo "Output: ${OUTDIR}"
echo "Started: $(date)"
echo ""

# Create output directories
mkdir -p "${OUTDIR}/subsampled"
mkdir -p "${OUTDIR}/full_db"
mkdir -p "${OUTDIR}/mtdna_db"

# Check databases exist
if [[ ! -f "${FULL_DB}/hash.k2d" ]]; then
    echo "ERROR: Full database not found at ${FULL_DB}"
    exit 1
fi

if [[ ! -f "${MTDNA_DB}/hash.k2d" ]]; then
    echo "ERROR: mtDNA database not found at ${MTDNA_DB}"
    echo "Please build it first with: sbatch build_mtdna_db.sh"
    exit 1
fi

# Parse samplesheet and process each sample
echo "=== Processing samples ==="
tail -n +2 "${SAMPLESHEET}" | while IFS=',' read -r sample fastq_1 fastq_2 sample_name species; do
    # Use sample_name if available, otherwise use sample
    name="${sample_name:-$sample}"

    echo ""
    echo "--- Processing: ${name} (${species}) ---"

    # Subsample reads
    echo "Subsampling to ${SUBSAMPLE_SIZE} reads..."
    subsampled="${OUTDIR}/subsampled/${name}.subsampled.fastq.gz"

    if [[ -n "$fastq_2" && -f "$fastq_2" ]]; then
        # Paired-end: subsample R1 only for simplicity
        run_seqtk seqtk sample -s42 "${fastq_1}" "${SUBSAMPLE_SIZE}" | gzip > "${subsampled}"
    else
        run_seqtk seqtk sample -s42 "${fastq_1}" "${SUBSAMPLE_SIZE}" | gzip > "${subsampled}"
    fi

    # Run Kraken2 with full database
    echo "Running Kraken2 with full database..."
    run_kraken2 kraken2 \
        --db "${FULL_DB}" \
        --threads ${SLURM_CPUS_PER_TASK:-8} \
        --report "${OUTDIR}/full_db/${name}.kraken2.report.txt" \
        --output "${OUTDIR}/full_db/${name}.kraken2.output.txt" \
        --gzip-compressed \
        --memory-mapping \
        "${subsampled}"

    # Run Kraken2 with mtDNA database
    echo "Running Kraken2 with mtDNA database..."
    run_kraken2 kraken2 \
        --db "${MTDNA_DB}" \
        --threads ${SLURM_CPUS_PER_TASK:-8} \
        --report "${OUTDIR}/mtdna_db/${name}.kraken2.report.txt" \
        --output "${OUTDIR}/mtdna_db/${name}.kraken2.output.txt" \
        --gzip-compressed \
        "${subsampled}"

    echo "Done with ${name}"
done

# Generate comparison report
echo ""
echo "=== Generating comparison report ==="

python3 << 'PYTHON_SCRIPT'
import os
import glob
import csv

outdir = os.environ.get('OUTDIR', 'kraken2_comparison')

def parse_kraken_report(report_file):
    """Parse Kraken2 report and extract key metrics."""
    percent_unclassified = 0.0
    top_species = "Unknown"
    top_species_percent = 0.0
    total_classified = 0

    if not os.path.exists(report_file):
        return None

    with open(report_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                percent = float(parts[0].strip())
                reads_clade = int(parts[1].strip())
                rank = parts[3].strip()
                taxon = parts[5].strip()

                if rank == 'U':
                    percent_unclassified = percent

                if rank == 'S' and percent > top_species_percent:
                    top_species_percent = percent
                    top_species = taxon

    percent_classified = 100.0 - percent_unclassified
    return {
        'percent_classified': percent_classified,
        'top_species': top_species,
        'top_species_percent': top_species_percent
    }

# Find all samples
full_reports = glob.glob(f"{outdir}/full_db/*.kraken2.report.txt")
samples = [os.path.basename(f).replace('.kraken2.report.txt', '') for f in full_reports]

# Parse results
results = []
for sample in sorted(samples):
    full_report = f"{outdir}/full_db/{sample}.kraken2.report.txt"
    mtdna_report = f"{outdir}/mtdna_db/{sample}.kraken2.report.txt"

    full_result = parse_kraken_report(full_report)
    mtdna_result = parse_kraken_report(mtdna_report)

    if full_result and mtdna_result:
        # Check if species match
        species_match = "YES" if full_result['top_species'] == mtdna_result['top_species'] else "NO"

        # Check if at least genus matches
        full_genus = full_result['top_species'].split()[0] if full_result['top_species'] != "Unknown" else ""
        mtdna_genus = mtdna_result['top_species'].split()[0] if mtdna_result['top_species'] != "Unknown" else ""
        genus_match = "YES" if full_genus == mtdna_genus and full_genus else "NO"

        results.append({
            'sample': sample,
            'full_pct_classified': full_result['percent_classified'],
            'full_top_species': full_result['top_species'],
            'full_top_pct': full_result['top_species_percent'],
            'mtdna_pct_classified': mtdna_result['percent_classified'],
            'mtdna_top_species': mtdna_result['top_species'],
            'mtdna_top_pct': mtdna_result['top_species_percent'],
            'species_match': species_match,
            'genus_match': genus_match
        })

# Print summary report
print("=" * 120)
print("KRAKEN2 DATABASE COMPARISON REPORT")
print("=" * 120)
print()

# Table header
header = f"{'Sample':<20} {'Full DB':<40} {'mtDNA DB':<40} {'Match':<10}"
print(header)
print(f"{'':20} {'% Class':<8} {'% Top':<8} {'Species':<24} {'% mtDNA':<8} {'% Top':<8} {'Species':<24} {'Sp/Gen':<10}")
print("-" * 120)

for r in results:
    full_species_short = r['full_top_species'][:22] + '..' if len(r['full_top_species']) > 24 else r['full_top_species']
    mtdna_species_short = r['mtdna_top_species'][:22] + '..' if len(r['mtdna_top_species']) > 24 else r['mtdna_top_species']

    print(f"{r['sample']:<20} {r['full_pct_classified']:>6.2f}%  {r['full_top_pct']:>6.2f}%  {full_species_short:<24} "
          f"{r['mtdna_pct_classified']:>6.2f}%  {r['mtdna_top_pct']:>6.2f}%  {mtdna_species_short:<24} "
          f"{r['species_match']}/{r['genus_match']}")

print()
print("=" * 120)
print("SUMMARY")
print("=" * 120)

total = len(results)
species_matches = sum(1 for r in results if r['species_match'] == 'YES')
genus_matches = sum(1 for r in results if r['genus_match'] == 'YES')
avg_mtdna = sum(r['mtdna_pct_classified'] for r in results) / total if total > 0 else 0

print(f"Total samples:        {total}")
print(f"Species match:        {species_matches}/{total} ({100*species_matches/total:.1f}%)")
print(f"Genus match:          {genus_matches}/{total} ({100*genus_matches/total:.1f}%)")
print(f"Avg % mtDNA reads:    {avg_mtdna:.2f}%")
print()

# Write CSV for further analysis
csv_file = f"{outdir}/comparison_results.csv"
with open(csv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)
print(f"Results saved to: {csv_file}")

# Write detailed report
report_file = f"{outdir}/comparison_report.txt"
with open(report_file, 'w') as f:
    f.write("KRAKEN2 DATABASE COMPARISON REPORT\n")
    f.write("=" * 80 + "\n\n")

    for r in results:
        f.write(f"Sample: {r['sample']}\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Full Database:\n")
        f.write(f"    % Classified:  {r['full_pct_classified']:.2f}%\n")
        f.write(f"    Top Species:   {r['full_top_species']} ({r['full_top_pct']:.2f}%)\n")
        f.write(f"  mtDNA Database:\n")
        f.write(f"    % mtDNA:       {r['mtdna_pct_classified']:.2f}%\n")
        f.write(f"    Top Species:   {r['mtdna_top_species']} ({r['mtdna_top_pct']:.2f}%)\n")
        f.write(f"  Species Match:   {r['species_match']}\n")
        f.write(f"  Genus Match:     {r['genus_match']}\n")
        f.write("\n")

print(f"Detailed report saved to: {report_file}")
PYTHON_SCRIPT

export OUTDIR="${OUTDIR}"

echo ""
echo "=== Comparison complete ==="
echo "Finished: $(date)"
echo ""
echo "Output files:"
echo "  - ${OUTDIR}/comparison_results.csv"
echo "  - ${OUTDIR}/comparison_report.txt"
echo "  - ${OUTDIR}/full_db/*.kraken2.report.txt"
echo "  - ${OUTDIR}/mtdna_db/*.kraken2.report.txt"
