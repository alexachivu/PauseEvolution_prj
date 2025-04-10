# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Load concatenated count data from intervals
load("/workdir/cgd24/NELF_dTAG_HS/interval_counts.RData")

# Concatenate counts from all intervals into one data frame
interval1_counts$region <- "interval1"
interval2_counts$region <- "interval2"
if (!is.null(interval3_counts)) {
    interval3_counts$region <- "interval3"
}

# Combine all intervals
combined_counts <- rbind(interval1_counts, interval2_counts, if (!is.null(interval3_counts)) interval3_counts else NULL)

# Remove gene_name and region columns to prepare for PCA
pca_data <- combined_counts[, !(names(combined_counts) %in% c("gene_name", "region"))]

# Filter out rows with zero counts across all samples
pca_data <- pca_data[rowSums(pca_data) > 0, ]

# Perform log transformation to stabilize variance
pca_data <- log2(pca_data + 1)

# Perform PCA
pca_result <- prcomp(t(pca_data), scale. = TRUE)

# Extract PC scores for plotting
pc_scores <- as.data.frame(pca_result$x)

# Extract sample information for coloring and point type
sample_info <- data.frame(sample_name = colnames(pca_data))
sample_info$NELF <- ifelse(grepl("NELFb", sample_info$sample_name), "NELFb", "NELFe")
sample_info$dTAG <- ifelse(grepl("NdT", sample_info$sample_name), "No dTAG", "dTAG")
sample_info$HS <- ifelse(grepl("NHS", sample_info$sample_name), "No Heat Stress", "Heat Stress")

# Merge PCA scores with sample information
pc_scores <- cbind(pc_scores, sample_info)

# Plot PCA results
p1 <- ggplot(pc_scores, aes(x = PC1, y = PC3, color = NELF, shape = interaction(dTAG, HS))) +
  geom_point(size = 3, aes(fill = interaction(dTAG, HS))) +
  theme_minimal() +
  labs(title = "PC1 vs PC3", x = "PC1", y = "PC3")

p2 <- ggplot(pc_scores, aes(x = PC2, y = PC3, color = NELF, shape = interaction(dTAG, HS))) +
  geom_point(size = 3, aes(fill = interaction(dTAG, HS))) +
  theme_minimal() +
  labs(title = "PC2 vs PC3", x = "PC2", y = "PC3")

p3 <- ggplot(pc_scores, aes(x = PC4, y = PC3, color = NELF, shape = interaction(dTAG, HS))) +
  geom_point(size = 3, aes(fill = interaction(dTAG, HS))) +
  theme_minimal() +
  labs(title = "PC4 vs PC3", x = "PC4", y = "PC3")

p4 <- ggplot(pc_scores, aes(x = PC5, y = PC3, color = NELF, shape = interaction(dTAG, HS))) +
  geom_point(size = 3, aes(fill = interaction(dTAG, HS))) +
  theme_minimal() +
  labs(title = "PC5 vs PC3", x = "PC5", y = "PC3")

p5 <- ggplot(pc_scores, aes(x = PC6, y = PC3, color = NELF, shape = interaction(dTAG, HS))) +
  geom_point(size = 3, aes(fill = interaction(dTAG, HS))) +
  theme_minimal() +
  labs(title = "PC6 vs PC3", x = "PC6", y = "PC3")

pc12 <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = NELF, shape = interaction(dTAG, HS))) +
	geom_point(size = 4, aes(fill = interaction(dTAG, HS))) +
	theme_minimal(base_size = 16) +  # Increases base font size
	labs(title = "PC1 vs PC2", x = "PC1", y = "PC2") +
	theme(
  		plot.title = element_text(size = 18, face = "bold"),
		axis.title = element_text(size = 16),
		axis.text = element_text(size = 14),
		legend.title = element_text(size = 14),
		legend.text = element_text(size = 12)
	)

# Save the plots to a PDF
pdf("/workdir/cgd24/NELF_dTAG_HS/PCA_plots.pdf", width = 25, height = 5)
grid.arrange(p1, p2, p3, p4, p5, ncol = 5)
dev.off()

pdf("~/transfer/NELF_dTAG_PC1_PC2.pdf")
grid.arrange(pc12)
dev.off()

cat("PCA plots saved to: /workdir/cgd24/NELF_dTAG_HS/PCA_plots.pdf\n")

