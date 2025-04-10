# Load necessary libraries
require(bigWig)
library(data.table)
library(GenomicRanges)

# Set file paths
bigwig_dir <- "/workdir/cgd24/NELF_dTAG_HS/aligned"
annotation_file <- "/local/storage/data/mm10/refGene.bed.gz"  # Corrected annotation file path

# Read gene annotations
annotations <- fread(cmd = paste("zcat", annotation_file), header = FALSE)
colnames(annotations) <- c("chr", "start", "end", "transcript_id", "score", "strand", "gene_name")

# Filter out poorly defined chromosomes (those with an underscore character)
annotations <- annotations[!grepl("_", annotations$chr), ]

# Filter transcripts shorter than 500 bp
annotations <- annotations[(annotations$end - annotations$start) >= 500, ]

# Convert annotations to GRanges
gene_ranges <- makeGRangesFromDataFrame(annotations,
                                        seqnames.field = "chr",
                                        start.field = "start",
                                        end.field = "end",
                                        strand.field = "strand",
                                        keep.extra.columns = TRUE)

# Vectorized function to split ranges into intervals
split_intervals_vectorized <- function(gr) {
    # Extract relevant fields
    tss <- ifelse(as.character(strand(gr)) == "+", start(gr), end(gr))
    tx_length <- width(gr)
    
    # Calculate interval boundaries for each range in a vectorized manner
    is_positive_strand <- as.character(strand(gr)) == "+"
    
    # Interval 1: 0-300 bp
    interval1_start <- ifelse(is_positive_strand, start(gr), pmax(end(gr) - 300, start(gr)))
    interval1_end <- ifelse(is_positive_strand, pmin(start(gr) + 300, end(gr)), end(gr))
    interval1 <- GRanges(seqnames(gr), IRanges(interval1_start, interval1_end), strand = strand(gr), gene_name = mcols(gr)$gene_name, score = mcols(gr)$score)
    
    # Interval 2: 300-3000 bp
    interval2_start <- ifelse(is_positive_strand, pmin(start(gr) + 300, end(gr)), pmax(end(gr) - 3000, start(gr)))
    interval2_end <- ifelse(is_positive_strand, pmin(start(gr) + 3000, end(gr)), pmax(end(gr) - 300, start(gr)))
    interval2 <- GRanges(seqnames(gr), IRanges(interval2_start, interval2_end), strand = strand(gr), gene_name = mcols(gr)$gene_name, score = mcols(gr)$score)
    
    # Interval 3: Remaining length after 3000 bp
    interval3_start <- ifelse(is_positive_strand, interval2_end + 1, start(gr))
    interval3_end <- ifelse(is_positive_strand, end(gr), interval2_start - 1)
    valid_interval3 <- ifelse(is_positive_strand, interval2_end < end(gr), interval2_start > start(gr))
    interval3 <- GRanges(seqnames(gr)[valid_interval3], IRanges(interval3_start[valid_interval3], interval3_end[valid_interval3]), strand = strand(gr)[valid_interval3], gene_name = mcols(gr)$gene_name[valid_interval3], score = mcols(gr)$score[valid_interval3])
    
    list(interval1 = interval1, interval2 = interval2, interval3 = interval3)
}

# Apply split_intervals_vectorized to all transcripts
split_list <- split_intervals_vectorized(gene_ranges)

# Extract intervals and convert to data frames
interval1 <- split_list$interval1
interval2 <- split_list$interval2
interval3 <- split_list$interval3

# Convert to BED6 format
interval1_bed6 <- as.data.frame(interval1)[, c("seqnames", "start", "end", "gene_name", "score", "strand")]
colnames(interval1_bed6) <- c("chr", "start", "end", "name", "score", "strand")

interval2_bed6 <- as.data.frame(interval2)[, c("seqnames", "start", "end", "gene_name", "score", "strand")]
colnames(interval2_bed6) <- c("chr", "start", "end", "name", "score", "strand")

if (length(interval3) > 0) {
    interval3_bed6 <- as.data.frame(interval3)[, c("seqnames", "start", "end", "gene_name", "score", "strand")]
    colnames(interval3_bed6) <- c("chr", "start", "end", "name", "score", "strand")
} else {
    interval3_bed6 <- NULL
}

# Remove ranges with 0 length (where start == end) from the interval bed6 data frames
interval1_bed6 <- interval1_bed6[interval1_bed6$start != interval1_bed6$end, ]
interval2_bed6 <- interval2_bed6[interval2_bed6$start != interval2_bed6$end, ]
if (!is.null(interval3_bed6)) {
    interval3_bed6 <- interval3_bed6[interval3_bed6$start != interval3_bed6$end, ]
}

# Remove exact duplicates (same start and end coordinates) from the interval bed6 data frames
interval1_bed6 <- interval1_bed6[!duplicated(interval1_bed6[, c("chr", "start", "end", "strand")]), ]
interval2_bed6 <- interval2_bed6[!duplicated(interval2_bed6[, c("chr", "start", "end", "strand")]), ]
if (!is.null(interval3_bed6)) {
    interval3_bed6 <- interval3_bed6[!duplicated(interval3_bed6[, c("chr", "start", "end", "strand")]), ]
}

# Save an image containing all three interval*_bed6 variables
save(interval1_bed6, interval2_bed6, interval3_bed6, file = "intervals_bed6.RData")

