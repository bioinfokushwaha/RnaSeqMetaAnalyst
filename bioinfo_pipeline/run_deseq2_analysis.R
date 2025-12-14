#!/usr/bin/env Rscript

# Complete DESeq2 Analysis with EnhancedVolcano Plots
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(EnhancedVolcano)
})

set.seed(0)

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_deseq2_analysis.R <count_file> <sampleinfo_file>")
}

count_file <- args[1]
sampleinfo_file <- args[2]

cat("\n==================================================\n")
cat("DESeq2 Differential Expression Analysis - COMPLETE\n")
cat("==================================================\n")
cat("Count file:", count_file, "\n")
cat("Sample info:", sampleinfo_file, "\n\n")

# Setup
filename <- tools::file_path_sans_ext(basename(count_file))
output_dir <- file.path("DEG", filename)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === 1. READ DATA ===
cat(">>> Step 1: Reading data...\n")
counts <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE, 
                   check.names = FALSE, row.names = 1)
cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

sampleinfo <- read.delim(sampleinfo_file, header = TRUE, sep = "\t", 
                         stringsAsFactors = FALSE)
colnames(sampleinfo) <- trimws(colnames(sampleinfo))
sampleinfo$SampleName <- trimws(sampleinfo$SampleName)
sampleinfo$CONDITION <- trimws(sampleinfo$CONDITION)

cat("  Sample info:", nrow(sampleinfo), "samples\n")
cat("  Conditions:", paste(unique(sampleinfo$CONDITION), collapse = ", "), "\n")

# === 2. MATCH SAMPLES ===
cat("\n>>> Step 2: Matching sample names...\n")
count_cols <- colnames(counts)
count_cols_clean <- gsub("\\.bam$|\\_Aligned.*$|\\.htseq$", "", count_cols)

sampleinfo$CleanName <- gsub("\\.bam$|\\_Aligned.*$|\\.fastq.*$", "", 
                             sampleinfo$SampleName)

matching <- intersect(count_cols_clean, sampleinfo$CleanName)
cat("  Matched samples:", length(matching), "\n")

if (length(matching) == 0) {
  stop("ERROR: No matching samples!")
}

# Reorder
count_idx <- match(sampleinfo$CleanName, count_cols_clean)
valid_idx <- !is.na(count_idx)
count_matrix <- as.matrix(counts[, count_idx[valid_idx]])
sampleinfo_matched <- sampleinfo[valid_idx, ]
count_matrix <- round(count_matrix)
mode(count_matrix) <- "integer"

cat("  Final dimensions:", nrow(count_matrix), "genes x", 
    ncol(count_matrix), "samples\n")
cat("  ℹ️  Decimal counts rounded to nearest integer\n") 

# === 3. FILTER LOW COUNTS ===
cat("\n>>> Step 3: Filtering low-count genes...\n")
genes_before <- nrow(count_matrix)
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix <- count_matrix[keep, ]
cat("  Kept", nrow(count_matrix), "of", genes_before, "genes\n")

# === 4. CREATE DESEQ2 OBJECT ===
cat("\n>>> Step 4: Creating DESeq2 object...\n")
colData <- data.frame(
  CONDITION = factor(sampleinfo_matched$CONDITION),
  row.names = colnames(count_matrix)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ CONDITION
)

# === 5. RUN DESEQ2 ===
cat("\n>>> Step 5: Running DESeq2...\n")
set.seed(0)
dds <- DESeq(dds)
cat("  ✓ DESeq2 completed\n")

# === 6. SAVE NORMALIZED COUNTS ===
cat("\n>>> Step 6: Saving normalized counts...\n")
normalized_counts <- counts(dds, normalized = TRUE)
normalized_df <- data.frame(
  Geneid = rownames(normalized_counts),
  normalized_counts,
  check.names = FALSE
)

norm_file <- file.path(output_dir, paste0(filename, "_DESeq2_Normalized_Counts.xlsx"))
write.xlsx(normalized_df, norm_file, rowNames = FALSE)
cat("  ✓ Saved:", norm_file, "\n")

# Create annotation with mean normalized counts
NorCount_anno <- data.frame(
  Trans = rownames(normalized_counts),
  NormalizedCount_Mean = rowMeans(normalized_counts),
  stringsAsFactors = FALSE
)

# === 7. HOUSEKEEPING GENES ANALYSIS ===
cat("\n>>> Step 7: Housekeeping genes analysis...\n")
hkg_file <- "/data/housekeeping_genes.xlsx"

if (file.exists(hkg_file)) {
  cat("  Found housekeeping genes file\n")
  
  tryCatch({
    hkg_data <- read.xlsx(hkg_file, sheet = 1)
    
    if (ncol(hkg_data) >= 1) {
      hkg_genes <- as.character(hkg_data[[1]])
      hkg_genes <- hkg_genes[!is.na(hkg_genes) & hkg_genes != ""]
      cat("  HKG list contains", length(hkg_genes), "genes\n")
      
      # Merge with normalized counts
      hkg_in_data <- hkg_genes[hkg_genes %in% normalized_df$Geneid]
      cat("  Found", length(hkg_in_data), "HKG in dataset\n")
      
      if (length(hkg_in_data) > 0) {
        hkg_normalized <- normalized_df[normalized_df$Geneid %in% hkg_in_data, ]
        
        # Filter for non-zero expression (all samples > 0)
        sample_cols <- 2:ncol(hkg_normalized)
        non_zero_mask <- rowSums(hkg_normalized[, sample_cols] > 0) == length(sample_cols)
        hkg_non_zero <- hkg_normalized[non_zero_mask, ]
        
        cat("  HKG with non-zero expression in all samples:", nrow(hkg_non_zero), "\n")
        
        if (nrow(hkg_non_zero) > 0) {
          hkg_out_file <- file.path(output_dir, 
                                     paste0(filename, "_DESeq2_HKG_NonZero_Normalized.xlsx"))
          write.xlsx(hkg_non_zero, hkg_out_file, rowNames = FALSE)
          cat("  ✓ Saved HKG non-zero expression file\n")
        } else {
          cat("  ⚠ No HKG with non-zero expression across all samples\n")
        }
      } else {
        cat("  ⚠ No HKG found in dataset\n")
      }
    }
  }, error = function(e) {
    cat("  ⚠ Error processing HKG file:", e$message, "\n")
  })
} else {
  cat("  ⚠ Housekeeping genes file not found at", hkg_file, "\n")
  cat("  Skipping HKG analysis\n")
}

# === 8. PCA ANALYSIS ===
cat("\n>>> Step 8: PCA analysis and plotting...\n")
set.seed(0)
vsd <- vst(dds, blind = TRUE)
pca_result <- prcomp(t(assay(vsd)))

# Calculate variance explained
pca_var <- pca_result$sdev^2
pca_var_pct <- round(100 * pca_var / sum(pca_var), 2)

# PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)
pca_scores$Condition <- colData$CONDITION

# Save PCA scores
pca_scores_file <- file.path(output_dir, paste0(filename, "_DESeq2_PCA_Scores.xlsx"))
write.xlsx(pca_scores, pca_scores_file, rowNames = FALSE)
cat("  ✓ Saved PCA scores\n")

# Create PCA plot (PC1 vs PC2)
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3, show.legend = FALSE) +
  labs(
    title = paste("PCA Plot -", filename, "- DESeq2"),
    x = paste0("PC1 (", pca_var_pct[1], "% variance)"),
    y = paste0("PC2 (", pca_var_pct[2], "% variance)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Save PCA plot
pca_plot_file <- file.path(output_dir, paste0(filename, "_DESeq2_PCA_Plot.png"))
ggsave(pca_plot_file, pca_plot, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved PCA plot (PC1 vs PC2)\n")

# === 9. DIFFERENTIAL EXPRESSION ANALYSIS WITH VOLCANO PLOTS ===
cat("\n>>> Step 9: Differential expression analysis and volcano plots...\n")
results_summary <- data.frame(
  Comparison = character(),
  UP_genes = numeric(),
  DOWN_genes = numeric(),
  Total_DEGs = numeric(),
  stringsAsFactors = FALSE
)

conditions <- levels(colData$CONDITION)
cat("  Conditions:", paste(conditions, collapse = ", "), "\n\n")

for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    comp_name <- paste0(cond1, "_vs_", cond2)
    cat("  Processing:", comp_name, "...\n")
    
    tryCatch({
      set.seed(0)
      res <- results(dds, contrast = c("CONDITION", cond1, cond2))
      res_df <- as.data.frame(res)
      res_df$Trans <- rownames(res_df)
      
      sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
      cat("    Significant DEGs:", nrow(sig), "\n")
      
      if (nrow(sig) > 0) {
        UP <- sig[sig$log2FoldChange >= 0, ]
        DOWN <- sig[sig$log2FoldChange < 0, ]
        
        cat("    UP:", nrow(UP), "| DOWN:", nrow(DOWN), "\n")
        
        # Format with normalized counts
        format_degs <- function(deg_df) {
          deg_formatted <- deg_df %>%
            select(Trans, logFC = log2FoldChange, logCPM = baseMean, 
                   LR = stat, PValue = pvalue, FDR = padj) %>%
            mutate(logCPM = log2(logCPM + 1))
          
          # Add normalized counts
          deg_anno <- left_join(deg_formatted, NorCount_anno, by = "Trans")
          
          # Filter for non-zero expression
          deg_non_zero <- deg_anno[deg_anno$NormalizedCount_Mean > 0, ]
          
          return(deg_non_zero)
        }
        
        if (nrow(UP) > 0) {
          UP_final <- format_degs(UP)
          if (nrow(UP_final) > 0) {
            up_file <- file.path(output_dir, 
                                paste0(filename, "_DESeq2_", comp_name, "_UP.xlsx"))
            write.xlsx(UP_final, up_file, rowNames = FALSE)
            cat("    ✓ Saved UP genes (non-zero):", nrow(UP_final), "\n")
          }
        }
        
        if (nrow(DOWN) > 0) {
          DOWN_final <- format_degs(DOWN)
          if (nrow(DOWN_final) > 0) {
            down_file <- file.path(output_dir, 
                                  paste0(filename, "_DESeq2_", comp_name, "_DOWN.xlsx"))
            write.xlsx(DOWN_final, down_file, rowNames = FALSE)
            cat("    ✓ Saved DOWN genes (non-zero):", nrow(DOWN_final), "\n")
          }
        }
        
        # === CREATE VOLCANO PLOT ===
        cat("    Creating EnhancedVolcano plot...\n")
        tryCatch({
          # Prepare data for volcano plot
          volcano_data <- res_df
          rownames(volcano_data) <- volcano_data$Trans
          
          # Create volcano plot
          volcano_plot <- EnhancedVolcano(
            volcano_data,
            lab = "",
            x = 'log2FoldChange',
            y = 'padj',
            title = paste(filename, "-", comp_name),
            subtitle = '',
            caption = "",
            legendLabels = c('NS', 
                           expression(log[2]~FC),
                           'adj p-value', 
                           expression(adj~p-value~~and~~log[2]~FC)),
            xlim = c(-8, 8),
            ylim = c(0, -log10(min(volcano_data$padj[volcano_data$padj > 0], na.rm = TRUE)) + 1),
            FCcutoff = 1.0,
            pCutoff = 0.05,
            pointSize = 3.0,
            labSize = 4.0,
            boxedLabels = FALSE,
            colAlpha = 4/5,
            legend = c('on', 'bottom'),
            legendPosition = 'bottom'
          )
          
          # Customize plot
          volcano_plot <- volcano_plot + ggplot2::theme(
            legend.text = ggplot2::element_text(size = 10),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = ggplot2::element_text(size = 12, face = "bold")
          )
          
          # Save volcano plot (high resolution)
          volcano_file <- file.path(output_dir, 
                                   paste0(filename, "_DESeq2_", comp_name, "_Volcano.png"))
          png(volcano_file, width = 3400, height = 3400, res = 600)
          print(volcano_plot)
          dev.off()
          
          cat("    ✓ Saved volcano plot:", volcano_file, "\n")
        }, error = function(e) {
          cat("    ✗ Error creating volcano plot:", e$message, "\n")
        })
        
        results_summary <- rbind(results_summary, data.frame(
          Comparison = comp_name,
          UP_genes = nrow(UP),
          DOWN_genes = nrow(DOWN),
          Total_DEGs = nrow(sig)
        ))
      } else {
        cat("    ⚠ No significant DEGs\n")
      }
    }, error = function(e) {
      cat("    ✗ ERROR:", e$message, "\n")
    })
  }
}

# === 10. SAVE SUMMARY ===
cat("\n>>> Step 10: Saving summary...\n")
summary_file <- file.path(output_dir, paste0(filename, "_DESeq2_DEG_Summary.xlsx"))
write.xlsx(results_summary, summary_file, rowNames = FALSE)

cat("\n==================================================\n")
cat("DESeq2 Analysis Complete!\n")
cat("==================================================\n")
cat("Output directory:", output_dir, "\n")
cat("Files generated:\n")
cat("  • Normalized counts\n")
cat("  • HKG non-zero expression (if applicable)\n")
cat("  • PCA scores and plot (PC1 vs PC2)\n")
cat("  • DEG lists (UP/DOWN, non-zero only)\n")
cat("  • Volcano plots (high resolution 3400x3400 @ 600 DPI)\n")
cat("  • Summary\n")
cat("\nTotal comparisons:", nrow(results_summary), "\n")
print(results_summary)
cat("\n")