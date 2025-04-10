# Load necessary libraries
library(Rsamtools)
library(data.table)
library(GenomicRanges)

# Set file paths
bam_dir <- "/workdir/cgd24/NELF_dTAG_HS/aligned/"

# List all bigWig files
bam_files <- list.files(bam_dir, pattern = "sorted.bam$", full.names = TRUE)

print(bam_files)

# A list to hold spike-in counts per sample
spikein_counts_list <- list()
mouse_counts_list <- list()

# Loop over BAM files and calculate number of reads for chroms containing 'dm3'
for (bam_file in bam_files) {
    print(bam_file)

    # Create a BamFile object and open it
    bam_obj <- BamFile(bam_file)
    open(bam_obj)
    
    # Get basic statistics for each reference sequence in the BAM
    idxstats <- idxstatsBam(bam_obj)
    
    # Close the BamFile
    close(bam_obj)
    
    # The columns in idxstats typically are: seqnames, seqlength, mapped, unmapped
    # Subset to chromosomes containing 'dm3' in their name
    dm3_rows <- grep("dm3", idxstats[,1])
    mm10_rows <- grep("dm3", idxstats[,1], invert=TRUE)
    
    if (length(dm3_rows) > 0) {
        # Sum all reads mapped to 'dm3' references
        total_spikein_reads <- sum(idxstats[dm3_rows, "mapped"])
    } else {
        total_spikein_reads <- 0
    }

    if(length(mm10_rows) > 0) {
	total_mouse_reads <- sum(idxstats[mm10_rows, "mapped"])
    } else {
	total_mouse_reads <- 0
    }
    
    # Store the spike-in counts, using the filename (sans extension) as a key
    sample_name <- sub("\\_QC_end.sort.sorted.bam$", "", basename(bam_file))
    spikein_counts_list[[sample_name]] <- total_spikein_reads
    mouse_counts_list[[sample_name]] <- total_mouse_reads
}

# Combine counts into a data frame
spikein_counts_df <- data.frame(
    sample = names(spikein_counts_list),
    spikein_reads = unlist(spikein_counts_list),
    mouse_reads = unlist(mouse_counts_list),
    row.names = NULL
)

spikein_counts_df

# Save the spike-in count table as RData (or CSV, etc.)
output_file_counts <- "/workdir/cgd24/NELF_dTAG_HS/spikein_counts.RData"
save(spikein_counts_df, file = output_file_counts)

cat("Spike-in count table saved to:", output_file_counts, "\n")


