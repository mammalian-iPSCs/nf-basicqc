#!/bin/bash
# Create samplesheet from CNAG XLS file
# Usage: ./create_samplesheet.sh <project.xls> <output.csv>
#
# The script auto-detects FASTQ paths from the XLS file location:
#   <project_dir>/<date_folder>/FASTQ/
#
# If multiple date folders exist, all are searched for each sample's FASTQs.
#
# Example:
#   ./create_samplesheet.sh /path/to/CGLZOO_01/CGLZOO_01.xls samplesheet.csv

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <project.xls> <output.csv>"
    echo ""
    echo "Arguments:"
    echo "  project.xls   - CNAG project Excel file"
    echo "  output.csv    - Output samplesheet path"
    echo ""
    echo "The FASTQ paths are auto-detected from: <project_dir>/<date_folder>/FASTQ/"
    echo "Multiple date folders are supported."
    echo ""
    echo "Example:"
    echo "  $0 /scratch/data_transfer/CGLZOO_01/CGLZOO_01.xls CGLZOO_01_samplesheet.csv"
    exit 1
fi

XLS_FILE="$1"
OUTPUT_CSV="$2"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if XLS file exists
if [ ! -f "$XLS_FILE" ]; then
    echo "Error: XLS file not found: $XLS_FILE"
    exit 1
fi

# Get project directory from XLS path
PROJECT_DIR="$(dirname "$XLS_FILE")"

# Find all date folders (directories with only digits) containing FASTQ subfolders
FASTQ_PATHS=()
for dir in "$PROJECT_DIR"/*/; do
    dirname=$(basename "$dir")
    if [[ "$dirname" =~ ^[0-9]+$ ]] && [ -d "${dir}FASTQ" ]; then
        FASTQ_PATHS+=("${dir}FASTQ")
    fi
done

if [ ${#FASTQ_PATHS[@]} -eq 0 ]; then
    echo "Error: Could not find any FASTQ directories in $PROJECT_DIR/<date>/FASTQ/"
    echo "Looking for numeric directories containing a FASTQ subfolder"
    exit 1
fi

echo "Found ${#FASTQ_PATHS[@]} FASTQ directory(ies):"
for path in "${FASTQ_PATHS[@]}"; do
    echo "  - $path"
done

# Join paths with comma for R script
FASTQ_PATHS_STR=$(IFS=','; echo "${FASTQ_PATHS[*]}")

# Load R module if on cluster
if command -v module &> /dev/null; then
    module load R 2>/dev/null || true
fi

# Run the R script
Rscript "${SCRIPT_DIR}/create_samplesheet.R" "$XLS_FILE" "$FASTQ_PATHS_STR" "$OUTPUT_CSV"
