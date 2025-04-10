# Load necessary libraries
library(DESeq2)

# Load concatenated count data from intervals
load("/workdir/cgd24/NELF_dTAG_HS/interval_counts.RData")

# Load spike-in counts
load("spikein_counts.RData")  # should create spikein_counts_df

# Focus on interval 3 for DESeq2 analysis
if (is.null(interval3_counts)) {
    stop("Interval 3 counts are not available.")
}

# Filter samples to include only those in the HS/NHS condition, with and without dTAG treatment
selected_samples <- grep("NELFb.(NdT|dT).*(HS|NHS)", colnames(interval3_counts))
filtered_counts <- interval3_counts[, selected_samples]

# Prepare unique gene names using a random number since start and end are not available
set.seed(42)  # For reproducibility
rownames(filtered_counts) <- paste0(interval3_counts[, "gene_name"], "_", sample(1e6, nrow(filtered_counts), replace = FALSE))

# Prepare colData for DESeq2, including sample conditions, cell lines, and dTAG treatment
sample_names <- colnames(filtered_counts)
condition <- ifelse(grepl("NHS", sample_names), "No_Heat_Shock", "Heat_Shock")
cell_line <- ifelse(grepl("NELFb", sample_names), "NELFb", "NELFe")
dtag_treatment <- ifelse(grepl("NdT", sample_names), "No_dTAG", "dTAG")
col_data <- data.frame(row.names = sample_names, 
                       condition = factor(condition), 
                       cell_line = factor(cell_line),
                       dtag_treatment = factor(dtag_treatment))

# Prepare DESeq2 dataset with condition, cell line, and dTAG treatment as fixed effects, including interaction term
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = col_data,
                              design = ~ dtag_treatment + condition)

dds$condition <- relevel(dds$condition, ref = "No_Heat_Shock")
dds$dtag_treatment <- relevel(dds$dtag_treatment, ref = "No_dTAG")

# --- Add spike-in size factors ---

# Ensure matching sample order
#stopifnot(all(spikein_counts_df$sample %in% colnames(dds))) ## Returns FALSE because we've already removed NELFe.
spikein_counts_df <- spikein_counts_df[match(colnames(dds), spikein_counts_df$sample), ]

# Compute inverse fraction as size factor, normalize to geometric mean 1
size_factors <- spikein_counts_df$spikein_reads
size_factors <- size_factors / mean(size_factors) #exp(mean(log(size_factors))) ## Alex used mean, not geometric mean.

# Set size factors manually
sizeFactors(dds) <- size_factors


# Run DESeq2 differential expression analysis. 
#dds <- DESeq(dds) ## DO NOT DO THIS since we are manually setting size factors.

# Run DESeq2, skipping size factor estimation
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

results_names <- resultsNames(dds)

# Extract interaction results for Heat Shock vs No Heat Shock modulated by dTAG treatment
res <- results(dds, name = "condition_Heat_Shock_vs_No_Heat_Shock")

# Sort results by adjusted p-value
res <- res[order(res$padj, na.last = NA), ]

# NOW: Take genes that are up-regulated in the condition condition_Heat_Shock_vs_No_Heat_Shock
upregulated_genes <- res[res$log2FoldChange > 0 & res$padj < 0.05, ]
upregulated_genes

downregulated_genes <- res[res$log2FoldChange < 0 & res$padj < 0.05, ]
downregulated_genes

## NOW Get adjusted counts for each condition, averaging the replicates.
# Extract model-adjusted counts for all samples
adjusted_counts <- counts(dds, normalized = TRUE)
adjusted_counts <- as.data.frame(adjusted_counts)

# Helper function to extract prefixes (first three parts of column names)
extract_prefix <- function(name) {
  sapply(strsplit(name, "\\."), function(x) paste(x[1:3], collapse = "."))
}

# Extract prefixes from column names
prefixes <- extract_prefix(colnames(adjusted_counts))
unique_prefixes <- unique(prefixes)  # Get unique prefixes

# Initialize an empty data frame to store results
averaged_counts <- data.frame(row.names = rownames(adjusted_counts))

# Loop over unique prefixes
for (prefix in unique_prefixes) {
  # Identify columns belonging to the current prefix
  matching_columns <- which(prefixes == prefix)
  
  # Calculate rowMeans for those columns
  averaged_counts[[prefix]] <- rowMeans(adjusted_counts[, matching_columns, drop = FALSE], na.rm = TRUE)
}

# Resulting data frame with averaged counts
head(averaged_counts)

## Now check out summaries for changes.
# Step 1: Filter columns containing "NELFb"
columns_nelfb <- grep("NELFb", colnames(averaged_counts), value = TRUE)
nelfb_counts <- averaged_counts[, columns_nelfb]

# Step 2: Filter rows where res$padj < 0.05 and res$log2FoldChange > 0.0
significant_genes_upregulated <- rownames(res[res$padj < 0.05 & res$log2FoldChange > 0.0, ])
filtered_counts <- nelfb_counts[rownames(nelfb_counts) %in% significant_genes_upregulated, ]

# Step 3: Subset rows with little/no difference between "NELFb.NdT.NHS" and "NELFb.dT.NHS"
# Retrieve the two specific columns
# Define a pseudocount
pc <- 0.01

ndt_nhs <- filtered_counts[, "NELFb.NdT.NHS"] + pc
dt_nhs <- filtered_counts[, "NELFb.dT.NHS"] + pc

# Identify genes with abs(log2 fold change) < 0.5 between these two columns
small_diff_genes <- abs(log2(ndt_nhs / dt_nhs)) < 0.5
filtered_counts <- filtered_counts[small_diff_genes, ]

# Step 4: Compute log2 fold change for specific column comparisons
ndt_hs <- filtered_counts[, "NELFb.NdT.HS"] + pc
dt_hs <- filtered_counts[, "NELFb.dT.HS"] + pc

ndt_nhs <- filtered_counts[, "NELFb.NdT.NHS"] + pc
dt_nhs <- filtered_counts[, "NELFb.dT.NHS"] + pc

# Compute log2 fold change
log2fc_ndt <- log2(ndt_hs / ndt_nhs)
log2fc_dt <- log2(dt_hs / dt_nhs)

# Generate summary statistics for the log2 fold changes
summary_ndt <- summary(log2fc_ndt)
summary_dt <- summary(log2fc_dt)

# Print results
cat("Summary of log2 fold change for NELFb.NdT:\n")
print(summary_ndt)

cat("\nSummary of log2 fold change for NELFb.dT:\n")
print(summary_dt)

cat("\nWilcox test for difference:\n")
wilcox.test(log2fc_ndt, log2fc_dt)

