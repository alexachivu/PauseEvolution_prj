# Load necessary libraries
require(bigWig)
library(data.table)
library(GenomicRanges)

# Load the intervals data from RData file
load("intervals_bed6.RData")

# Set file paths
bigwig_dir <- "/workdir/cgd24/NELF_dTAG_HS/aligned/My_output-01_08_2025/"

# List all bigWig files
bigwig_files <- list.files(bigwig_dir, pattern = ".bw$", full.names = TRUE)

# Split files by strand
bigwig_plus <- bigwig_files[grepl("_plus.bw$", bigwig_files)]
bigwig_minus <- bigwig_files[grepl("_minus.bw$", bigwig_files)]

# Ensure strand consistency
stopifnot(length(bigwig_plus) == length(bigwig_minus))

# Load bigWig files
bw_files <- list()
for (i in seq_along(bigwig_plus)) {
    sample_name <- sub("_QC_end.sort_plus.bw$", "", basename(bigwig_plus[i]))
    bw_files[[paste0(sample_name, "_pl")]] <- load.bigWig(bigwig_plus[i])
    bw_files[[paste0(sample_name, "_mn")]] <- load.bigWig(bigwig_minus[i])
}

# Get counts for each sample using bed6.region.bpQuery.bigWig
count_matrix_interval1 <- list()
count_matrix_interval2 <- list()
count_matrix_interval3 <- list()

for (sample_name in names(bw_files)) {
    if (grepl("_pl$", sample_name)) {
        sample_root <- sub("_pl$", "", sample_name)
        plus_bw <- bw_files[[sample_name]]
        minus_bw <- bw_files[[paste0(sample_root, "_mn")]]

        # Get read counts for interval 1, combining plus and minus strands
        counts_interval1 <- bed6.region.bpQuery.bigWig(plus_bw, minus_bw, interval1_bed6, abs.value = TRUE)
        count_matrix_interval1[[sample_root]] <- counts_interval1

        # Get read counts for interval 2, combining plus and minus strands
        counts_interval2 <- bed6.region.bpQuery.bigWig(plus_bw, minus_bw, interval2_bed6, abs.value = TRUE)
        count_matrix_interval2[[sample_root]] <- counts_interval2

        # Get read counts for interval 3, combining plus and minus strands (if interval3 exists)
        if (!is.null(interval3_bed6)) {
            counts_interval3 <- bed6.region.bpQuery.bigWig(plus_bw, minus_bw, interval3_bed6, abs.value = TRUE)
            count_matrix_interval3[[sample_root]] <- counts_interval3
        }
    }
}

# Combine all counts into data frames
interval1_counts <- do.call(cbind, count_matrix_interval1)
interval1_counts <- data.frame(gene_name = interval1_bed6$name, interval1_counts)

interval2_counts <- do.call(cbind, count_matrix_interval2)
interval2_counts <- data.frame(gene_name = interval2_bed6$name, interval2_counts)

if (!is.null(interval3_bed6)) {
    interval3_counts <- do.call(cbind, count_matrix_interval3)
    interval3_counts <- data.frame(gene_name = interval3_bed6$name, interval3_counts)
} else {
    interval3_counts <- NULL
}

# Save the count matrices as a single RData file
output_file_counts <- "/workdir/cgd24/NELF_dTAG_HS/interval_counts.RData"
save(interval1_counts, interval2_counts, interval3_counts, file = output_file_counts)

cat("Count tables saved to:", output_file_counts, "\n")

# Clean up: unload bigWig files
for (sample_name in names(bw_files)) {
    unload.bigWig(bw_files[[sample_name]])
}

