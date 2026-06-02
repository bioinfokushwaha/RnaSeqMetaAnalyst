#!/usr/bin/env Rscript
# ===================================================================
# edgeR Differential Expression Analysis - Complete Version
# ===================================================================
# Features: Project directory support, robust sample matching,
# EnhancedVolcano plots, HKG analysis with PCA, comprehensive logging
# ===================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(EnhancedVolcano)
})

set.seed(0)

# -------------------------------------------------------------------
# PARSE ARGUMENTS
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_edgeR_analysis.R <count_file> <sampleinfo_file> [--project-dir <path>] [--log2fc x --fdr y --pvalue z]")
}

count_file <- args[1]
sampleinfo_file <- args[2]

# Default values
LOG2FC_THRESHOLD <- 1.0
FDR_THRESHOLD <- 0.05
PVALUE_THRESHOLD <- 0.05
PROJECT_DIR <- "/data"  # FIXED: Always use /data (the volume-mounted directory)

# Parse optional arguments
i <- 3
while (i <= length(args)) {
  if (args[i] == "--project-dir") {
    PROJECT_DIR <- args[i+1]
    i <- i + 2
  } else if (args[i] == "--log2fc") {
    LOG2FC_THRESHOLD <- as.numeric(args[i+1])
    i <- i + 2
  } else if (args[i] == "--fdr") {
    FDR_THRESHOLD <- as.numeric(args[i+1])
    i <- i + 2
  } else if (args[i] == "--pvalue") {
    PVALUE_THRESHOLD <- as.numeric(args[i+1])
    i <- i + 2
  } else {
    i <- i + 1
  }
}

# -------------------------------------------------------------------
# SETUP OUTPUT DIRECTORY
# -------------------------------------------------------------------
file_prefix <- sub("_counts_clean\\.txt$", "", basename(count_file))

# Create organized directory structure: /data/results/DEG/<pipeline>/edgeR/
pipeline_subdir <- file.path(PROJECT_DIR, "results", "DEG", file_prefix)
dir.create(pipeline_subdir, showWarnings = FALSE, recursive = TRUE)

output_dir <- file.path(pipeline_subdir, "edgeR")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Verify directory was created and is writable
if (!dir.exists(output_dir)) {
  stop(paste("CRITICAL ERROR: Failed to create output directory:", output_dir))
}
if (!file.access(output_dir, mode = 2) == 0) {
  stop(paste("CRITICAL ERROR: Output directory not writable:", output_dir))
}

# Setup logging
log_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_Analysis.log"))
log_conn <- file(log_file, "w")
cat_log <- function(...) { 
  msg <- paste(...)
  cat(msg)
  cat(msg, file=log_conn)
}

cat_log("=======================================================\n")
cat_log("edgeR Differential Expression Analysis\n")
cat_log("=======================================================\n")
cat_log("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat_log("Project Directory: ", PROJECT_DIR, "\n")
cat_log("Output Directory: ", output_dir, "\n")
cat_log("ABSOLUTE PATH: ", normalizePath(output_dir, mustWork = FALSE), "\n")
cat_log("Directory exists: ", dir.exists(output_dir), "\n")
cat_log("Directory writable: ", file.access(output_dir, mode = 2) == 0, "\n")
cat_log("Count File: ", count_file, "\n")
cat_log("Sample Info: ", sampleinfo_file, "\n")
cat_log("Thresholds:\n")
cat_log("  - log2FC: ", LOG2FC_THRESHOLD, "\n")
cat_log("  - FDR: ", FDR_THRESHOLD, "\n")
cat_log("  - p-value: ", PVALUE_THRESHOLD, "\n")
cat_log("=======================================================\n\n")

# -------------------------------------------------------------------
# STEP 1: READ DATA
# -------------------------------------------------------------------
cat_log(">>> Step 1: Reading data...\n")
counts <- read.csv(count_file, header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE, row.names = 1)
cat_log("  Count matrix: ", nrow(counts), " genes x ", ncol(counts), " samples\n")

sampleinfo <- read.delim(sampleinfo_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

# Flexible column name detection
sample_col <- ifelse("SampleName" %in% colnames(sampleinfo), "SampleName",
              ifelse("Sample" %in% colnames(sampleinfo), "Sample", colnames(sampleinfo)[1]))
cond_col <- ifelse("CONDITION" %in% colnames(sampleinfo), "CONDITION",
            ifelse("Condition" %in% colnames(sampleinfo), "Condition", colnames(sampleinfo)[2]))

# Clean whitespace
colnames(sampleinfo) <- trimws(colnames(sampleinfo))
sampleinfo[[sample_col]] <- trimws(sampleinfo[[sample_col]])
sampleinfo[[cond_col]] <- trimws(sampleinfo[[cond_col]])

cat_log("  Sample info: ", nrow(sampleinfo), " samples\n")
cat_log("  Conditions: ", paste(unique(sampleinfo[[cond_col]]), collapse=", "), "\n\n")

# -------------------------------------------------------------------
# STEP 2: ROBUST SAMPLE NAME MATCHING
# -------------------------------------------------------------------
cat_log(">>> Step 2: Matching sample names...\n")
count_cols <- colnames(counts)

# Comprehensive sample name cleaning function
clean_sample_name <- function(name) {
  name <- trimws(name)
  name <- gsub("\\.bam$", "", name)
  name <- gsub("_Aligned\\..*", "", name)
  name <- gsub("\\.htseq$", "", name)
  name <- gsub("\\.htseq\\.txt$", "", name)
  name <- gsub("\\.genes\\.results$", "", name)
  name <- gsub("\\.fastq.*$", "", name)
  name <- gsub("_trimmed$", "", name)
  name <- gsub("_R1$|_R2$", "", name)
  name <- gsub("_sorted$", "", name)
  name <- trimws(name)
  return(name)
}

count_cols_clean <- sapply(count_cols, clean_sample_name)
sampleinfo$CleanName <- sapply(sampleinfo[[sample_col]], clean_sample_name)

cat_log("  Count matrix samples:\n")
cat_log("    Original: ", paste(count_cols[1:min(3, length(count_cols))], collapse=", "), "...\n")
cat_log("    Cleaned:  ", paste(count_cols_clean[1:min(3, length(count_cols_clean))], collapse=", "), "...\n\n")

cat_log("  Sample info names:\n")
cat_log("    Original: ", paste(sampleinfo[[sample_col]][1:min(3, nrow(sampleinfo))], collapse=", "), "...\n")
cat_log("    Cleaned:  ", paste(sampleinfo$CleanName[1:min(3, nrow(sampleinfo))], collapse=", "), "...\n\n")

# Match samples
matching <- intersect(count_cols_clean, sampleinfo$CleanName)
cat_log("  Matched samples: ", length(matching), "\n")

if (length(matching) == 0) {
  cat_log("\n  ERROR: No matching samples found!\n")
  cat_log("  Count matrix has: ", paste(count_cols_clean, collapse=", "), "\n")
  cat_log("  Sample info has: ", paste(sampleinfo$CleanName, collapse=", "), "\n")
  close(log_conn)
  stop("Sample name matching failed - check file names and sampleinfo.txt")
}

# Reorder and filter
count_idx <- match(sampleinfo$CleanName, count_cols_clean)
valid_idx <- !is.na(count_idx)

if (sum(valid_idx) == 0) {
  close(log_conn)
  stop("No valid samples after matching!")
}

count_matrix <- as.matrix(counts[, count_idx[valid_idx]])
sampleinfo_matched <- sampleinfo[valid_idx, ]
cat_log("  Valid matched samples: ", sum(valid_idx), "\n\n")

# -------------------------------------------------------------------
# STEP 3: HANDLE DECIMAL VALUES (RSEM COMPATIBILITY)
# -------------------------------------------------------------------
cat_log(">>> Step 3: Processing count values...\n")

has_decimals <- any(count_matrix != round(count_matrix), na.rm = TRUE)

if (has_decimals) {
  cat_log("  Detected decimal counts (likely RSEM output)\n")
  cat_log("  Rounding to nearest integer for edgeR\n")
  count_matrix <- round(count_matrix)
} else {
  cat_log("  ✓ Counts are already integers\n")
}

mode(count_matrix) <- "integer"
cat_log("  ✓ Converted to integer mode\n\n")

# -------------------------------------------------------------------
# STEP 4: CREATE EDGER OBJECT & FILTER
# -------------------------------------------------------------------
cat_log(">>> Step 4: Creating DGEList and filtering...\n")

dge <- DGEList(
  counts = count_matrix,
  group = factor(sampleinfo_matched[[cond_col]])
)

genes_before <- nrow(dge)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
cat_log("  Kept ", nrow(dge), " of ", genes_before, " genes\n")
cat_log("  Removed ", genes_before - nrow(dge), " low-count genes\n\n")

# -------------------------------------------------------------------
# STEP 5: NORMALIZE & ESTIMATE DISPERSION
# -------------------------------------------------------------------
cat_log(">>> Step 5: Normalization and dispersion estimation...\n")

dge <- calcNormFactors(dge)
cat_log("  ✓ Calculated normalization factors\n")

design <- model.matrix(~ group, data = dge$samples)
colnames(design) <- gsub("group", "", colnames(design))

dge <- estimateDisp(dge, design)
cat_log("  ✓ Estimated dispersion\n\n")

# -------------------------------------------------------------------
# STEP 6: SAVE NORMALIZED COUNTS
# -------------------------------------------------------------------
cat_log(">>> Step 6: Saving normalized counts...\n")

normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
normalized_df <- data.frame(
  Gene = rownames(normalized_counts),
  normalized_counts,
  check.names = FALSE
)

# Save Excel
norm_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_Normalized_Counts.xlsx"))
write.xlsx(normalized_df, norm_file, rowNames = FALSE)
cat_log("  ✓ Saved: ", basename(norm_file), "\n")
cat_log("     Full path: ", normalizePath(norm_file, mustWork = FALSE), "\n")

# Save CSV
norm_csv <- file.path(output_dir, paste0(file_prefix, "_edgeR_normalized.csv"))
write.csv(normalized_df, norm_csv, row.names = FALSE)
cat_log("  ✓ Saved: ", basename(norm_csv), "\n\n")

# Create mean normalized counts annotation
NorCount_anno <- data.frame(
  Gene = rownames(normalized_counts),
  NormalizedCount_Mean = rowMeans(normalized_counts),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------
# STEP 7: HOUSEKEEPING GENES ANALYSIS
# -------------------------------------------------------------------
cat_log(">>> Step 7: Housekeeping genes analysis...\n")
hkg_file <- file.path(PROJECT_DIR, "housekeeping_genes.xlsx")

if (file.exists(hkg_file)) {
  cat_log("  Found: ", hkg_file, "\n")
  
  tryCatch({
    hkg_data <- read.xlsx(hkg_file, sheet = 1)
    
    if (ncol(hkg_data) >= 1) {
      hkg_genes <- as.character(hkg_data[[1]])
      hkg_genes <- hkg_genes[!is.na(hkg_genes) & hkg_genes != ""]
      cat_log("  HKG list: ", length(hkg_genes), " genes\n")
      
      hkg_in_data <- hkg_genes[hkg_genes %in% normalized_df$Gene]
      cat_log("  Found in dataset: ", length(hkg_in_data), " HKG\n")
      
      if (length(hkg_in_data) > 0) {
        hkg_normalized <- normalized_df[normalized_df$Gene %in% hkg_in_data, ]
        
        # Filter non-zero in all samples
        sample_cols <- 2:ncol(hkg_normalized)
        non_zero_mask <- rowSums(hkg_normalized[, sample_cols] > 0) == length(sample_cols)
        hkg_non_zero <- hkg_normalized[non_zero_mask, ]
        
        cat_log("  Non-zero in all samples: ", nrow(hkg_non_zero), " HKG\n")
        
        if (nrow(hkg_non_zero) > 0) {
          # Save HKG normalized counts
          hkg_out_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_HKG_NonZero_Normalized.xlsx"))
          write.xlsx(hkg_non_zero, hkg_out_file, rowNames = FALSE)
          cat_log("  ✓ Saved: ", basename(hkg_out_file), "\n")
          
          # HKG PCA
          cat_log("  Generating HKG PCA...\n")
          log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
          hkg_genes_in_data <- rownames(log_cpm)[rownames(log_cpm) %in% hkg_in_data]
          
          if (length(hkg_genes_in_data) > 1) {
            hkg_log_cpm <- log_cpm[hkg_genes_in_data, ]
            hkg_pca <- prcomp(t(hkg_log_cpm))
            
            hkg_pca_var <- hkg_pca$sdev^2
            hkg_pca_var_pct <- round(100 * hkg_pca_var / sum(hkg_pca_var), 2)
            
            hkg_pca_df <- data.frame(
              Sample = rownames(hkg_pca$x),
              Condition = dge$samples$group,
              PC1 = hkg_pca$x[,1],
              PC2 = hkg_pca$x[,2]
            )
            
            # Plot
            hkg_pca_plot <- ggplot(hkg_pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
              geom_point(size = 4, alpha = 0.8) +
              geom_text(vjust = -0.8, hjust = 0.5, size = 3, show.legend = FALSE) +
              labs(
                title = paste("HKG PCA -", file_prefix),
                subtitle = paste("Based on", length(hkg_genes_in_data), "housekeeping genes"),
                x = paste0("PC1 (", hkg_pca_var_pct[1], "%)"),
                y = paste0("PC2 (", hkg_pca_var_pct[2], "%)")
              ) +
              theme_bw() +
              theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                plot.subtitle = element_text(hjust = 0.5, size = 10),
                legend.position = "right"
              )
            
            hkg_pca_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_HKG_PCA_Plot.png"))
            ggsave(hkg_pca_file, hkg_pca_plot, width = 10, height = 8, dpi = 300)
            cat_log("  ✓ Saved: ", basename(hkg_pca_file), "\n")
            
            # Save scores
            hkg_scores_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_HKG_PCA_Scores.xlsx"))
            write.xlsx(hkg_pca_df, hkg_scores_file, rowNames = FALSE)
            cat_log("  ✓ Saved: ", basename(hkg_scores_file), "\n")
          }
        }
      }
    }
  }, error = function(e) {
    cat_log("  ⚠ Error: ", e$message, "\n")
  })
} else {
  cat_log("  ⚠ Not found: ", hkg_file, "\n")
  cat_log("  Skipping HKG analysis\n")
}

cat_log("\n")

# -------------------------------------------------------------------
# STEP 8: PCA ANALYSIS (ALL GENES)
# -------------------------------------------------------------------
cat_log(">>> Step 8: PCA analysis...\n")
set.seed(0)
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
pca_result <- prcomp(t(log_cpm))

pca_var <- pca_result$sdev^2
pca_var_pct <- round(100 * pca_var / sum(pca_var), 2)

pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)
pca_scores$Condition <- dge$samples$group

# Save scores
pca_scores_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_PCA_Scores.xlsx"))
write.xlsx(pca_scores, pca_scores_file, rowNames = FALSE)
cat_log("  ✓ Saved PCA scores\n")

# Plot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3, show.legend = FALSE) +
  labs(
    title = paste("PCA -", file_prefix, "- edgeR"),
    x = paste0("PC1 (", pca_var_pct[1], "%)"),
    y = paste0("PC2 (", pca_var_pct[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

pca_plot_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_PCA_Plot.png"))
ggsave(pca_plot_file, pca_plot, width = 10, height = 8, dpi = 300)
cat_log("  ✓ Saved PCA plot\n\n")

# -------------------------------------------------------------------
# STEP 9: DIFFERENTIAL EXPRESSION WITH VOLCANO PLOTS
# -------------------------------------------------------------------
cat_log(">>> Step 9: Differential expression analysis...\n")

# Fit model
fit <- glmQLFit(dge, design)

results_summary <- data.frame(
  Comparison = character(),
  UP_genes = numeric(),
  DOWN_genes = numeric(),
  Total_DEGs = numeric(),
  stringsAsFactors = FALSE
)

conditions <- levels(dge$samples$group)
cat_log("  Conditions: ", paste(conditions, collapse=", "), "\n\n")

for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    cond1 <- conditions[i]
    cond2 <- conditions[j]
    comp_name <- paste0(cond1, "_vs_", cond2)
    cat_log("  Processing: ", comp_name, "\n")
    
    tryCatch({
      # Find contrast column
      contrast_col <- which(colnames(design) == cond2)
      
      if (length(contrast_col) == 0) {
        cat_log("    ⚠ Cannot find contrast column\n")
        next
      }
      
      set.seed(0)
      qlf <- glmQLFTest(fit, coef = contrast_col)
      res <- topTags(qlf, n = Inf, sort.by = "none")$table
      res$Gene <- rownames(res)
      
      # Apply thresholds
      sig <- res[res$FDR < FDR_THRESHOLD & abs(res$logFC) >= LOG2FC_THRESHOLD, ]
      
      cat_log("    Significant: ", nrow(sig), " DEGs\n")
      
      if (nrow(sig) > 0) {
        UP <- sig[sig$logFC > 0, ]
        DOWN <- sig[sig$logFC < 0, ]
        
        cat_log("    UP: ", nrow(UP), " | DOWN: ", nrow(DOWN), "\n")
        
        # Format function
        # FIX: glmQLFTest returns 'F' statistic column, NOT 'LR'.
        # Selecting 'LR' caused a silent error in tryCatch, skipping all saves.
        # Also removed NormalizedCount_Mean > 0 filter that silently dropped all DEGs.
        format_degs <- function(deg_df) {
          stat_col <- intersect(c("F", "LR"), colnames(deg_df))
          keep_cols <- c("Gene", "logFC", "logCPM", stat_col, "PValue", "FDR")
          keep_cols <- keep_cols[keep_cols %in% colnames(deg_df)]
          deg_formatted <- deg_df[, keep_cols, drop = FALSE]
          deg_anno <- left_join(deg_formatted, NorCount_anno, by = "Gene")
          return(deg_anno)
        }
        
        # Save UP
        if (nrow(UP) > 0) {
          UP_final <- format_degs(UP)
          up_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_", comp_name, "_UP.xlsx"))
          write.xlsx(UP_final, up_file, rowNames = FALSE)
          cat_log("    ✓ UP saved: ", nrow(UP_final), " genes\n")
          cat_log("       Saved to: ", normalizePath(up_file, mustWork = FALSE), "\n")
        }
        
        # Save DOWN
        if (nrow(DOWN) > 0) {
          DOWN_final <- format_degs(DOWN)
          down_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_", comp_name, "_DOWN.xlsx"))
          write.xlsx(DOWN_final, down_file, rowNames = FALSE)
          cat_log("    ✓ DOWN saved: ", nrow(DOWN_final), " genes\n")
          cat_log("       Saved to: ", normalizePath(down_file, mustWork = FALSE), "\n")
        }

        # Save all results
        all_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_", comp_name, "_All_Results.xlsx"))
        write.xlsx(res, all_file, rowNames = FALSE)
        
        # EnhancedVolcano plot
        cat_log("    Creating volcano plot...\n")
        tryCatch({
          volcano_data <- res
          rownames(volcano_data) <- volcano_data$Gene
          
          volcano_plot <- EnhancedVolcano(
            volcano_data,
            lab = "",
            x = 'logFC',
            y = 'FDR',
            title = paste(file_prefix, "-", comp_name),
            subtitle = paste0("edgeR | FC≥", LOG2FC_THRESHOLD, " | FDR<", FDR_THRESHOLD),
            caption = "",
            legendLabels = c('NS', expression(log[2]~FC), 'FDR', expression(FDR~~and~~log[2]~FC)),
            xlim = c(-8, 8),
            ylim = c(0, -log10(min(volcano_data$FDR[volcano_data$FDR > 0], na.rm = TRUE)) + 1),
            FCcutoff = LOG2FC_THRESHOLD,
            pCutoff = FDR_THRESHOLD,
            pointSize = 3.0,
            labSize = 4.0,
            boxedLabels = FALSE,
            colAlpha = 4/5,
            legendPosition = 'bottom',
            legendLabSize = 10,
            legendIconSize = 4.0
          )
          
          volcano_plot <- volcano_plot + theme(
            legend.text = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            axis.title = element_text(size = 12, face = "bold")
          )
          
          volcano_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_", comp_name, "_Volcano.png"))
          png(volcano_file, width = 3400, height = 3400, res = 600)
          print(volcano_plot)
          dev.off()
          
          cat_log("    ✓ Volcano: ", basename(volcano_file), "\n")
        }, error = function(e) {
          cat_log("    ✗ Volcano error: ", e$message, "\n")
        })
        
        results_summary <- rbind(results_summary, data.frame(
          Comparison = comp_name,
          UP_genes = nrow(UP),
          DOWN_genes = nrow(DOWN),
          Total_DEGs = nrow(sig)
        ))
      }
    }, error = function(e) {
      cat_log("    ✗ ERROR: ", e$message, "\n")
    })
  }
}

cat_log("\n")

# -------------------------------------------------------------------
# STEP 10: SAVE SUMMARY
# -------------------------------------------------------------------
cat_log(">>> Step 10: Saving summary...\n")
summary_file <- file.path(output_dir, paste0(file_prefix, "_edgeR_Summary.xlsx"))
write.xlsx(results_summary, summary_file, rowNames = FALSE)
cat_log("  ✓ Saved: ", basename(summary_file), "\n")

# -------------------------------------------------------------------
# COMPLETION
# -------------------------------------------------------------------
cat_log("\n=======================================================\n")
cat_log("edgeR Analysis Complete!\n")
cat_log("=======================================================\n")
cat_log("Output: ", output_dir, "\n")
cat_log("ABSOLUTE OUTPUT PATH: ", normalizePath(output_dir, mustWork = FALSE), "\n")
cat_log("Total comparisons: ", nrow(results_summary), "\n")
cat_log("Total DEGs: ", sum(results_summary$Total_DEGs), "\n")

if (nrow(results_summary) > 0) {
  cat_log("\nSummary:\n")
  cat_log(paste(capture.output(print(results_summary, row.names = FALSE)), collapse = "\n"), "\n")
}

cat_log("=======================================================\n")

close(log_conn)
cat("✓ edgeR completed\n")
cat("Results written to: ", normalizePath(output_dir, mustWork = FALSE), "\n")
