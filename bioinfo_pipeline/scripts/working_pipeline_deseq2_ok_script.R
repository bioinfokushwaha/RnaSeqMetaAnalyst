# Load libraries
library(DESeq2)
library(dplyr)
library(openxlsx)
set.seed(123)
# Read data
counts <- read.delim("gene_count_matrix.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sampleinfo <- read.delim("sampleinfo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Prepare count matrix
count_matrix <- counts[, -1]
rownames(count_matrix) <- counts$gene_id
count_matrix <- count_matrix[, sampleinfo$SampleName]  # Match order

# Convert to integer matrix (required for DESeq2)
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer"

# Create DESeq2 dataset
colData <- data.frame(
  SampleName = sampleinfo$SampleName,
  CONDITION = factor(sampleinfo$CONDITION),
  row.names = sampleinfo$SampleName
)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ CONDITION)

# Run DESeq2 normalization and differential expression
dds <- DESeq(dds)

# Get normalized counts for annotation
normalized_counts <- counts(dds, normalized = TRUE)
NorCount_anno <- data.frame(
  Trans = rownames(normalized_counts),
  NormalizedCount_Mean = rowMeans(normalized_counts),
  stringsAsFactors = FALSE
)

# Create output directory
dir.create("DEG_Results_DESeq2", showWarnings = FALSE)

# Initialize summary tracker
results_summary <- data.frame(
  Comparison = character(),
  UP_genes = numeric(),
  DOWN_genes = numeric(),
  Total_DEGs = numeric(),
  stringsAsFactors = FALSE
)

# List of unique conditions
condition <- unique(sampleinfo$CONDITION)

# Loop over condition pairs
for (i in 1:length(condition)) {
  for (j in (i + 1):length(condition)) {
    if (j <= length(condition)) {
      
      cond1 <- condition[i]
      cond2 <- condition[j]
      comp_name <- paste0(cond1, "_vs_", cond2)
      print(paste("Processing:", comp_name))
      
      tryCatch({
        # Get results for this comparison
        res <- results(dds, contrast = c("CONDITION", cond1, cond2))
        
        # Convert to data frame and add gene names
        res_df <- as.data.frame(res)
        res_df$Trans <- rownames(res_df)
        
        # Filter significant DEGs (padj < 0.05, equivalent to FDR < 0.05)
        sig <- res_df[which(!is.na(res_df$padj) & res_df$padj < 0.05), ]
        print(paste("Total significant DEGs:", nrow(sig)))
        
        if (nrow(sig) > 0) {
          # UP and DOWN (log2FoldChange >= 0 for UP, < 0 for DOWN)
          UP <- subset(sig, log2FoldChange >= 0)
          DOWN <- subset(sig, log2FoldChange < 0)
          
          print(paste("UP:", nrow(UP), "DOWN:", nrow(DOWN)))
          
          # Rename columns to match edgeR output format
          UP_formatted <- UP %>%
            select(Trans, logFC = log2FoldChange, logCPM = baseMean, 
                   LR = stat, PValue = pvalue, FDR = padj) %>%
            mutate(logCPM = log2(logCPM + 1))  # Convert to log scale
          
          DOWN_formatted <- DOWN %>%
            select(Trans, logFC = log2FoldChange, logCPM = baseMean, 
                   LR = stat, PValue = pvalue, FDR = padj) %>%
            mutate(logCPM = log2(logCPM + 1))  # Convert to log scale
          
          # Annotate with normalized counts
          UP_anno <- left_join(UP_formatted, NorCount_anno, by = "Trans")
          DOWN_anno <- left_join(DOWN_formatted, NorCount_anno, by = "Trans")
          
          # Save UP
          up_file <- "DEG_Results_DESeq2/Treatment_Rgene_UP_DESeq2.xlsx"
          wb_up <- if (file.exists(up_file)) loadWorkbook(up_file) else createWorkbook()
          addWorksheet(wb_up, comp_name)
          writeData(wb_up, comp_name, UP_anno, rowNames = FALSE)
          saveWorkbook(wb_up, up_file, overwrite = TRUE)
          
          # Save DOWN
          down_file <- "DEG_Results_DESeq2/Treatment_Rgene_DOWN_DESeq2.xlsx"
          wb_down <- if (file.exists(down_file)) loadWorkbook(down_file) else createWorkbook()
          addWorksheet(wb_down, comp_name)
          writeData(wb_down, comp_name, DOWN_anno, rowNames = FALSE)
          saveWorkbook(wb_down, down_file, overwrite = TRUE)
          
          # Track summary
          results_summary <- rbind(results_summary, data.frame(
            Comparison = comp_name,
            UP_genes = nrow(UP),
            DOWN_genes = nrow(DOWN),
            Total_DEGs = nrow(sig)
          ))
        } else {
          print("No significant DEGs")
          results_summary <- rbind(results_summary, data.frame(
            Comparison = comp_name,
            UP_genes = 0,
            DOWN_genes = 0,
            Total_DEGs = 0
          ))
        }
        
      }, error = function(e) {
        print(paste("Error in", comp_name, ":", e$message))
      })
      
      print("---")
    }
  }
}

# Save summary
write.xlsx(results_summary, "DEG_Results_DESeq2/DEG_Analysis_Summary_DESeq2.xlsx", rowNames = FALSE)

# Print summary
print("=== DESEQ2 ANALYSIS COMPLETE ===")
print(results_summary)
