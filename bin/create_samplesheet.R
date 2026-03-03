#!/usr/bin/env Rscript
# Script to create a samplesheet from CNAG XLS file
# Usage: Rscript create_samplesheet.R <project.xls> <fastq_paths> <output.csv>
#
# fastq_paths can be comma-separated if multiple directories exist

library(readxl)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript create_samplesheet.R <project.xls> <fastq_paths> <output.csv>")
}

xls_path <- args[1]
fastq_paths_str <- args[2]
output_csv <- args[3]

# Parse comma-separated FASTQ paths
fastq_paths <- unlist(strsplit(fastq_paths_str, ","))
cat("Searching for FASTQs in", length(fastq_paths), "directory(ies)\n")

# Read the Excel file (skip first 2 rows as in original script)
xls_data <- read_excel(xls_path, skip = 2)

# Function to find FASTQ file across multiple directories
find_fastq <- function(sample_id, suffix, paths) {
  filename <- paste0(sample_id, suffix)
  for (path in paths) {
    full_path <- file.path(path, filename)
    if (file.exists(full_path)) {
      return(full_path)
    }
  }
  # Return NA if not found
  warning(paste("FASTQ not found:", filename))
  return(NA)
}

# Create samplesheet dataframe
samplesheet <- data.frame(
  sample = xls_data$FLI,
  fastq_1 = sapply(xls_data$FLI, function(x) find_fastq(x, "_1.fastq.gz", fastq_paths)),
  fastq_2 = sapply(xls_data$FLI, function(x) find_fastq(x, "_2.fastq.gz", fastq_paths)),
  sample_name = xls_data$`SAMPLE BARCODE`,
  species = xls_data$SPECIES,
  stringsAsFactors = FALSE
)

# Check for missing files
missing <- sum(is.na(samplesheet$fastq_1))
if (missing > 0) {
  cat("Warning:", missing, "samples have missing FASTQ files\n")
}

# Write output CSV
write.csv(samplesheet, output_csv, row.names = FALSE, quote = FALSE)

cat("Samplesheet created:", output_csv, "\n")
cat("Number of samples:", nrow(samplesheet), "\n")
