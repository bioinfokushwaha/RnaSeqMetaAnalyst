#!/usr/bin/env Rscript

# Complete edgeR Analysis with EnhancedVolcano Plots
suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(EnhancedVolcano)
})

set.seed(0)

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_edgeR_analysis.R <count_file> <sampleinfo_file>")
}

count_file <- args[1]
sampleinfo_file <- args[2]

cat("\n==================================================\n")
cat("edgeR Differential Expression Analysis - COMPLETE\n")
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

cat("  Final dimensions:", nrow(count_matrix), "genes x", 
    ncol(count_matrix), "samples\n")

# === 3. CREATE DGELIST & FILTER ===
cat("\n>>> Step 3: Creating DGEList and filtering...\n")
dge <- DGEList(counts = count_matrix, group = sampleinfo_matched$CONDITION)

genes_before <- nrow(dge)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat("  Kept", nrow(dge), "of", genes_before, "genes\n")

# === 4. NORMALIZE ===
cat("\n>>> Step 4: Normalizing...\n")
dge <- calcNormFactors(dge)
cat("  ✓ Normalization completed\n")

# === 5. SAVE NORMALIZED COUNTS (CPM) ===
cat("\n>>> Step 5: Saving normalized counts...\n")
cpm_normalized <- cpm(dge, normalized.lib.sizes = TRUE)
normalized_df <- data.frame(
  Geneid = rownames(cpm_normalized),
  cpm_normalized,
  check.names = FALSE
)

norm_file <- file.path(output_dir, paste0(filename, "_edgeR_Normalized_CPM.xlsx"))
write.xlsx(normalized_df, norm_file, rowNames = FALSE)
cat("  ✓ Saved:", norm_file, "\n")

# Create annotation with mean normalized CPM
NorCount_anno <- data.frame(
  Trans = rownames(cpm_normalized),
  NormalizedCount_Mean = rowMeans(cpm_normalized),
  stringsAsFactors = FALSE
)

# === 6. HOUSEKEEPING GENES ANALYSIS ===
cat("\n>>> Step 6: Housekeeping genes analysis...\n")
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
                                     paste0(filename, "_edgeR_HKG_NonZero_Normalized.xlsx"))
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

# === 7. PCA ANALYSIS ===
cat("\n>>> Step 7: PCA analysis and plotting...\n")
set.seed(0)
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
pca_result <- prcomp(t(logCPM))

# Calculate variance explained
pca_var <- pca_result$sdev^2
pca_var_pct <- round(100 * pca_var / sum(pca_var), 2)

# PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)
pca_scores$Condition <- sampleinfo_matched$CONDITION

# Save PCA scores
pca_scores_file <- file.path(output_dir, paste0(filename, "_edgeR_PCA_Scores.xlsx"))
write.xlsx(pca_scores, pca_scores_file, rowNames = FALSE)
cat("  ✓ Saved PCA scores\n")

# Create PCA plot (PC1 vs PC2)
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3, show.legend = FALSE) +
  labs(
    title = paste("PCA Plot -", filename, "- edgeR"),
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
pca_plot_file <- file.path(output_dir, paste0(filename, "_edgeR_PCA_Plot.png"))
ggsave(pca_plot_file, pca_plot, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved PCA plot (PC1 vs PC2)\n")

# === 8. DIFFERENTIAL EXPRESSION ANALYSIS WITH VOLCANO PLOTS ===
cat("\n>>> Step 8: Differential expression analysis and volcano plots...\n")
results_summary <- data.frame(
  Comparison = character(),
  UP_genes = numeric(),
  DOWN_genes = numeric(),
  Total_DEGs = numeric(),
  stringsAsFactors = FALSE
)

conditions <- unique(sampleinfo_matched$CONDITION)
cat("  Conditions:", paste(conditions, collapse = ", "), "\n\n")

for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    comp_name <- paste0(cond1, "_vs_", cond2)
    cat("  Processing:", comp_name, "...\n")
    
    tryCatch({
      keep_samples <- sampleinfo_matched$CONDITION %in% c(cond1, cond2)
      dge_sub <- dge[, keep_samples]
      dge_sub$samples$group <- droplevels(factor(sampleinfo_matched$CONDITION[keep_samples]))
      
      design <- model.matrix(~0 + dge_sub$samples$group)
      colnames(design) <- levels(dge_sub$samples$group)
      
      set.seed(0)
      dge_sub <- estimateDisp(dge_sub, design)
      fit <- glmFit(dge_sub, design)
      
      contrast <- makeContrasts(
        contrasts = paste0("`", cond1, "` - `", cond2, "`"),
        levels = design
      )
      
      set.seed(0)
      lrt <- glmLRT(fit, contrast = contrast)
      
      res <- topTags(lrt, n = Inf)$table
      res$Trans <- rownames(res)
      
      sig <- res[res$FDR < 0.05, ]
      cat("    Significant DEGs:", nrow(sig), "\n")
      
      if (nrow(sig) > 0) {
        UP <- sig[sig$logFC >= 0, ]
        DOWN <- sig[sig$logFC < 0, ]
        
        cat("    UP:", nrow(UP), "| DOWN:", nrow(DOWN), "\n")
        
        # Format with normalized CPM
        format_degs <- function(deg_df) {
          deg_formatted <- deg_df %>%
            select(Trans, logFC, logCPM, LR, PValue, FDR)
          
          # Add normalized CPM mean
          deg_anno <- left_join(deg_formatted, NorCount_anno, by = "Trans")
          
          # Filter for non-zero expression
          deg_non_zero <- deg_anno[deg_anno$NormalizedCount_Mean > 0, ]
          
          return(deg_non_zero)
        }
        
        if (nrow(UP) > 0) {
          UP_final <- format_degs(UP)
          if (nrow(UP_final) > 0) {
            up_file <- file.path(output_dir, 
                                paste0(filename, "_edgeR_", comp_name, "_UP.xlsx"))
            write.xlsx(UP_final, up_file, rowNames = FALSE)
            cat("    ✓ Saved UP genes (non-zero):", nrow(UP_final), "\n")
          }
        }
        
        if (nrow(DOWN) > 0) {
          DOWN_final <- format_degs(DOWN)
          if (nrow(DOWN_final) > 0) {
            down_file <- file.path(output_dir, 
                                  paste0(filename, "_edgeR_", comp_name, "_DOWN.xlsx"))
            write.xlsx(DOWN_final, down_file, rowNames = FALSE)
            cat("    ✓ Saved DOWN genes (non-zero):", nrow(DOWN_final), "\n")
          }
        }
        
        # === CREATE VOLCANO PLOT ===
        cat("    Creating EnhancedVolcano plot...\n")
        tryCatch({
          # Prepare data for volcano plot - convert to data frame format for EnhancedVolcano
          volcano_data <- res
          rownames(volcano_data) <- volcano_data$Trans
          
          # Create volcano plot
          volcano_plot <- EnhancedVolcano(
            volcano_data,
            lab = "",
            x = 'logFC',
            y = 'FDR',
            title = paste(filename, "-", comp_name),
            subtitle = '',
            caption = "",
            legendLabels = c('NS', 
                           expression(log[2]~FC),
                           'adj p-value', 
                           expression(adj~p-value~~and~~log[2]~FC)),
            xlim = c(-8, 8),
            ylim = c(0, -log10(min(volcano_data$FDR[volcano_data$FDR > 0], na.rm = TRUE)) + 1),
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
                                   paste0(filename, "_edgeR_", comp_name, "_Volcano.png"))
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

# === 9. SAVE SUMMARY ===
cat("\n>>> Step 9: Saving summary...\n")
summary_file <- file.path(output_dir, paste0(filename, "_edgeR_DEG_Summary.xlsx"))
write.xlsx(results_summary, summary_file, rowNames = FALSE)

cat("\n==================================================\n")
cat("edgeR Analysis Complete!\n")
cat("==================================================\n")
cat("Output directory:", output_dir, "\n")
cat("Files generated:\n")
cat("  • Normalized CPM\n")
cat("  • HKG non-zero expression (if applicable)\n")
cat("  • PCA scores and plot (PC1 vs PC2)\n")
cat("  • DEG lists (UP/DOWN, non-zero only)\n")
cat("  • Volcano plots (high resolution 3400x3400 @ 600 DPI)\n")
cat("  • Summary\n")
cat("\nTotal comparisons:", nrow(results_summary), "\n")
print(results_summary)
cat("\n")