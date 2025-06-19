#!/usr/bin/env Rscript

# Load libraries
library(DESeq2)
library(dplyr)
library(openxlsx)
library(tools)

# Command-line args: 1 = count file, 2 = sampleinfo file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_deseq2_analysis.R <count_file> <sampleinfo_file>")
}

count_file <- args[1]
sampleinfo_file <- args[2]

# Extract base name for output directory (without extension)
prefix <- file_path_sans_ext(basename(count_file))
output_dir <- file.path("DEG", prefix)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
counts <- read.delim(count_file, header = TRUE, stringsAsFactors = FALSE)
sampleinfo <- read.delim(sampleinfo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Prepare count matrix
count_matrix <- counts[, -1]
rownames(count_matrix) <- counts$Geneid
count_matrix <- count_matrix[, sampleinfo$SampleName]
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer"

# Build DESeq2 dataset
colData <- data.frame(
  SampleName = sampleinfo$SampleName,
  CONDITION = factor(sampleinfo$CONDITION),
  row.names = sampleinfo$SampleName
)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~CONDITION)
dds <- DESeq(dds)

# Normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
NorCount_anno <- data.frame(
  Trans = rownames(normalized_counts),
  NormalizedCount_Mean = rowMeans(normalized_counts),
  stringsAsFactors = FALSE
)

# Initialize summary
results_summary <- data.frame(Comparison = character(), UP_genes = numeric(),
                              DOWN_genes = numeric(), Total_DEGs = numeric(), stringsAsFactors = FALSE)

# Differential Expression Analysis
condition <- unique(sampleinfo$CONDITION)
for (i in 1:(length(condition)-1)) {
  for (j in (i+1):length(condition)) {
    cond1 <- condition[i]
    cond2 <- condition[j]
    comp_name <- paste0(cond1, "_vs_", cond2)
    message("Processing: ", comp_name)

    tryCatch({
      res <- results(dds, contrast = c("CONDITION", cond1, cond2))
      res_df <- as.data.frame(res)
      res_df$Trans <- rownames(res_df)
      sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

      if (nrow(sig) > 0) {
        UP <- subset(sig, log2FoldChange >= 0)
        DOWN <- subset(sig, log2FoldChange < 0)

        UP_fmt <- UP %>%
          select(Trans, logFC = log2FoldChange, logCPM = baseMean, LR = stat, PValue = pvalue, FDR = padj) %>%
          mutate(logCPM = log2(logCPM + 1))
        DOWN_fmt <- DOWN %>%
          select(Trans, logFC = log2FoldChange, logCPM = baseMean, LR = stat, PValue = pvalue, FDR = padj) %>%
          mutate(logCPM = log2(logCPM + 1))

        UP_anno <- left_join(UP_fmt, NorCount_anno, by = "Trans")
        DOWN_anno <- left_join(DOWN_fmt, NorCount_anno, by = "Trans")

        # Save results
        up_file <- file.path(output_dir, paste0(prefix, "_UP_DESeq2.xlsx"))
        wb_up <- if (file.exists(up_file)) loadWorkbook(up_file) else createWorkbook()
        addWorksheet(wb_up, comp_name)
        writeData(wb_up, comp_name, UP_anno, rowNames = FALSE)
        saveWorkbook(wb_up, up_file, overwrite = TRUE)

        down_file <- file.path(output_dir, paste0(prefix, "_DOWN_DESeq2.xlsx"))
        wb_down <- if (file.exists(down_file)) loadWorkbook(down_file) else createWorkbook()
        addWorksheet(wb_down, comp_name)
        writeData(wb_down, comp_name, DOWN_anno, rowNames = FALSE)
        saveWorkbook(wb_down, down_file, overwrite = TRUE)

        results_summary <- rbind(results_summary, data.frame(
          Comparison = comp_name,
          UP_genes = nrow(UP),
          DOWN_genes = nrow(DOWN),
          Total_DEGs = nrow(sig)
        ))
      } else {
        results_summary <- rbind(results_summary, data.frame(
          Comparison = comp_name, UP_genes = 0, DOWN_genes = 0, Total_DEGs = 0
        ))
      }
    }, error = function(e) {
      message("Error in ", comp_name, ": ", e$message)
    })
  }
}

# Save summary
write.xlsx(results_summary, file.path(output_dir, paste0(prefix, "_Summary_DESeq2.xlsx")), rowNames = FALSE)
message("=== DESEQ2 ANALYSIS COMPLETE ===")
print(results_summary)
