# =========================================================================
# 1. Load the libraries
# =========================================================================
library(DESeq2)
library(ggplot2)

# =========================================================================
# 2. Dataset downloaded and working directory changed to where the file is located
# =========================================================================

# Read the extracted CSV file 
excel_data <- read.csv("GSE326385_FFA_gene.csv", header = TRUE, check.names = FALSE)
counts_matrix <- as.data.frame(excel_data)

# Assign Column A containing Gene Symbols as the row index registry
rownames(counts_matrix) <- counts_matrix[, 1]

# Strip away Column A so only pure numeric sample count values remain
counts_matrix <- counts_matrix[, -1]

# Coerce counts into pure matrix integers to meet Bioconductor constraints
counts_matrix <- as.matrix(sapply(counts_matrix, as.integer))

# =========================================================================
# 3. Explicit Sample Mapping Metadata (18 Total Samples)
# =========================================================================
sample_names <- colnames(counts_matrix)

# Map the 9 Control samples followed directly by the 9 FFA-Treated samples 
conditions <- c(rep("Control", 9), rep("FFA_Treated", 9))

# Construct a metadata dataframe for statistical grouping pairing
metadata <- data.frame(sample_id = sample_names,condition = factor(conditions), 
                            levels = c("Control", "FFA_Treated"))
rownames(metadata) <- metadata$sample_id

# =========================================================================
# 4. Construct DESeq2 Object and Calculate Statistics
# =========================================================================
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,colData = metadata,
                              design = ~ condition)

# Filter low-expression background noise genes (minimum 10 cumulative reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run the core empirical differential expression calculation pipeline
dds <- DESeq(dds)

# Extract standard statistical metrics contrasting Treated vs Control Baseline
res <- results(dds, contrast = c("condition", "FFA_Treated", "Control"))
res_df <- as.data.frame(res)

# Export final tabular spreadsheet file for assignment delivery
write.csv(res_df, file = "GSE326385_DGE_Results.csv")
print("Task 1 Complete! Output table successfully saved as 'GSE326385_DGE_Results.csv'.")

# =========================================================================
# 5. Extract Specific Target Gene Values 
# =========================================================================
# Querying EGR2 expression statistics for verification
if("EGR2" %in% rownames(res_df)) {
  cat("\n--- EGR2 Target Gene Metrics discovered in your data ---\n")
  print(res_df["EGR2", ])
} else {
  cat("\nEGR2 not found. Showing the single highest differentially expressed gene instead:\n")
  print(head(res_df[order(res_df$padj), ], 1))
}

# =========================================================================
# 6. Visualization using Volcano Plot
# =========================================================================
# Categorize entries using strict significance criteria (p-adj < 0.05 & |Log2FC| > 1)
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Upregulated"
res_df$Significance[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Downregulated"

# Generate the plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: FFA Treatment vs Control (GSE326385)",
       subtitle = "Significance Thresholds: p-adjusted < 0.05 | |Log2FC| > 1",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))