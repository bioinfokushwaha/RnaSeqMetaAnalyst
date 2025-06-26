#!/usr/bin/env Rscript

# Load libraries
library(DESeq2)
library(dplyr)
library(openxlsx)

# Read command line argument
args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
sampleinfo_file <- "DEG/sampleinfo.txt"  # Adjust path if needed

# Parse filename without extension
filename <- tools::file_path_sans_ext(basename(count_file))
output_dir <- file.path("DEG", filename)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read data
counts <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE)
sampleinfo <- read.delim(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Prepare count matrix
count_matrix <- counts[, -1]
rownames(count_matrix) <- counts[, 1]
count_matrix <- count_matrix[, sampleinfo$SampleName]  # Match order

# Convert to integer matrix
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer"

# Create DESeq2 object
colData <- data.frame(
  SampleName = sampleinfo$SampleName,
  CONDITION = factor(sampleinfo$CONDITION),
  row.names = sampleinfo$SampleName
)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~CONDITION)
dds <- DESeq(dds)

# Normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
NorCount_anno <- data.frame(Trans = rownames(normalized_counts),
                            NormalizedCount_Mean = rowMeans(normalized_counts),
                            stringsAsFactors = FALSE)

# Summary holder
results_summary <- data.frame(Comparison = character(), UP_genes = numeric(), DOWN_genes = numeric(), Total_DEGs = numeric(), stringsAsFactors = FALSE)

condition <- unique(sampleinfo$CONDITION)

for (i in 1:(length(condition) - 1)) {
  for (j in (i + 1):length(condition)) {
    cond1 <- condition[i]
    cond2 <- condition[j]
    comp_name <- paste0(cond1, "_vs_", cond2)
    print(paste("Processing:", comp_name))

    tryCatch({
      res <- results(dds, contrast = c("CONDITION", cond1, cond2))
      res_df <- as.data.frame(res)
      res_df$Trans <- rownames(res_df)
      sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

      if (nrow(sig) > 0) {
        UP <- subset(sig, log2FoldChange >= 0)
        DOWN <- subset(sig, log2FoldChange < 0)

        # Format
        UP_formatted <- UP %>%
          select(Trans, logFC = log2FoldChange, logCPM = baseMean, LR = stat, PValue = pvalue, FDR = padj) %>%
          mutate(logCPM = log2(logCPM + 1))
        DOWN_formatted <- DOWN %>%
          select(Trans, logFC = log2FoldChange, logCPM = baseMean, LR = stat, PValue = pvalue, FDR = padj) %>%
          mutate(logCPM = log2(logCPM + 1))

        # Annotate
        UP_anno <- left_join(UP_formatted, NorCount_anno, by = "Trans")
        DOWN_anno <- left_join(DOWN_formatted, NorCount_anno, by = "Trans")

        # Save to Excel
        up_file <- file.path(output_dir, "UP_DESeq2.xlsx")
        down_file <- file.path(output_dir, "DOWN_DESeq2.xlsx")

        wb_up <- createWorkbook()
        addWorksheet(wb_up, comp_name)
        writeData(wb_up, comp_name, UP_anno)
        saveWorkbook(wb_up, up_file, overwrite = TRUE)

        wb_down <- createWorkbook()
        addWorksheet(wb_down, comp_name)
        writeData(wb_down, comp_name, DOWN_anno)
        saveWorkbook(wb_down, down_file, overwrite = TRUE)

        results_summary <- rbind(results_summary, data.frame(
          Comparison = comp_name,
          UP_genes = nrow(UP),
          DOWN_genes = nrow(DOWN),
          Total_DEGs = nrow(sig)
        ))
      }
    }, error = function(e) {
      message("Error in ", comp_name, ": ", e$message)
    })
  }
}

# Save summary
write.xlsx(results_summary, file.path(output_dir, "DEG_Analysis_Summary_DESeq2.xlsx"), rowNames = FALSE)
print("=== DESEQ2 ANALYSIS COMPLETE ===")
print(results_summary)
