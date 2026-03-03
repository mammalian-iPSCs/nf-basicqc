# Session Notes

## Current Status (2026-03-03)

### Working Features
- **Pipeline runs successfully** with all modules (FastQC, FastQ Screen, Kraken2, Sex Determination, SortMeRNA, RiboDetector, MultiQC)
- **Summary table output** (`results/summary/qc_summary.tsv`) - consolidated TSV with all QC metrics
- **Kraken2 plot without unclassified** - modified reports exclude unclassified reads, percentages recalculated
- **SortMeRNA rRNA quantification** - % rRNA in subsampled reads, MultiQC general stats integration
- **RiboDetector rRNA quantification** - % rRNA via deep learning, MultiQC general stats integration
- **rRNA Kraken2 (SILVA SSU)** - species ID from rRNA reads using SILVA 138 SSU database (new)

### rRNA Kraken2 Feature (added 2026-03-03)
- SortMeRNA conditionally saves rRNA reads when `--rrna_kraken2_db` is provided
- KRAKEN2_RRNA runs on those reads against a SILVA SSU Kraken2 database
- Results (top species, % of classified) appear as `rrna_top_species` / `rrna_pct_top_species` in `qc_summary.tsv`
- rRNA reads published to `results/sortmerna/` as `{fli}_rrna_1.fastq.gz` / `{fli}_rrna_2.fastq.gz`
- Kraken2 reports published to `results/kraken2_rrna/` as `{fli}_rrna.kraken2.report.txt`
- SILVA 138 SSU pre-built database: https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz (~112 MB)

## Key Files

### Pipeline Configuration
- `nextflow.config` - Container versions, resource settings, params
- `submit_pipeline.sh` - Production run script
- `test/submit_tests.sh` - Test scripts (use `--rrna_kraken2` to test new feature)

### Modules
- `modules/fastqc.nf` - FastQC
- `modules/fastq_screen.nf` - FastQ Screen
- `modules/kraken2.nf` - Kraken2 (used for both mtDNA and rRNA classification via alias)
- `modules/sex_determination.nf` - Sex determination + SUMMARIZE_SEX
- `modules/sortmerna.nf` - SortMeRNA (index + classification; now optionally saves rRNA reads)
- `modules/ribodetector.nf` - RiboDetector
- `modules/summarize_kraken2.nf` - Kraken2 summary + modified reports for MultiQC
- `modules/summarize_results.nf` - Consolidated QC summary table (qc_summary.tsv)
- `modules/prepare_multiqc_config.nf` - Generates MultiQC config YAML
- `modules/multiqc.nf` - MultiQC

### Custom Content Files Generated
- `kraken2_top_species_mqc.txt` - General stats columns (plot_type: 'generalstats')
- `*_classified.kraken2.report.txt` - Modified Kraken reports without unclassified
- `sex_determination_mqc.txt` - Sex determination results

## Databases
- Kraken2 mtDNA: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_mtdna`
- Sex markers: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/sex_markers/all_sex_markers.fasta`
- SortMeRNA rRNA FASTAs: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/rRNA_indices`
- SILVA SSU Kraken2 (to download): `https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz`
  - Suggested path: `/scratch_isilon/groups/compgen/data/Illumina_CryoZoo/genomes/databases/kraken/k2_silva_ssu`
  - Update `SILVA_DB` in `test/submit_tests.sh` once downloaded

### Callithrix species in mtDNA database
| Species | Taxon ID | K-mers |
|---------|----------|--------|
| C. aurita | 57375 | 2,240 |
| C. geoffroyi | 52231 | 853 |
| C. kuhlii | 867363 | 572 |
| C. jacchus | 9483 | 476 |
| C. penicillata | 57378 | 350 |

~39% of Callithrix k-mers are shared at genus level.

## Commands

### Run tests
```bash
sbatch test/submit_tests.sh --full             # Full pipeline
sbatch test/submit_tests.sh --sortmerna        # SortMeRNA rRNA
sbatch test/submit_tests.sh --rrna_kraken2     # rRNA Kraken2 (SILVA SSU) — needs SILVA_DB set
sbatch test/submit_tests.sh --kraken-fresh     # Kraken2 only, no resume
```

### Production run
```bash
sbatch submit_pipeline.sh <samplesheet.csv> <output_dir> [project_name] [application]
```

### Download SILVA SSU database
```bash
wget https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz
mkdir -p k2_silva_ssu && tar -xzf 16S_Silva138_20200326.tgz -C k2_silva_ssu/
```

## Recent Commits
```
Add rRNA-based species ID via Kraken2 on SortMeRNA output (SILVA SSU)
e4a179b Add SortMeRNA and RiboDetector % rRNA to consolidated summary table
b53d36c Add RiboDetector module and per-database SortMeRNA test
828ba3b Add SortMeRNA % rRNA to MultiQC general stats
194d827 Add standalone RiboDetector test script
```
