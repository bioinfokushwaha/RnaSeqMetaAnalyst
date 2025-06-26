#!/usr/bin/env Rscript

# Load libraries
library(edgeR)
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
count_matrix <- count_matrix[, sampleinfo$SampleName]  # Match sample order

# Create DGEList
dge <- DGEList(counts = count_matrix, group = sampleinfo$CONDITION)
dge <- calcNormFactors(dge)

# Normalized CPM for annotation
cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE)
NorCount_anno <- data.frame(
  Trans = rownames(cpm_norm),
  NormalizedCount_Mean = rowMeans(cpm_norm),
  stringsAsFactors = FALSE
)

# Summary holder
results_summary <- data.frame(Comparison = character(), UP_genes = numeric(), DOWN_genes = numeric(), Total_DEGs = numeric(), stringsAsFactors = FALSE)

condition <- unique(sampleinfo$CONDITION)

for (i in 1:(length(condition) - 1)) {
  for (j in (i + 1):length(condition)) {
    cond1 <- condition[i]
    cond2 <- condition[j]
    comp_name <- paste0(cond1, "_vs_", cond2)
    message("Processing: ", comp_name)

    tryCatch({
      keep_samples <- sampleinfo$CONDITION %in% c(cond1, cond2)
      dge_sub <- dge[, keep_samples]
      dge_sub$samples$group <- droplevels(factor(sampleinfo$CONDITION[keep_samples]))

      design <- model.matrix(~0 + dge_sub$samples$group)
      colnames(design) <- levels(dge_sub$samples$group)

      dge_sub <- estimateDisp(dge_sub, design)
      fit <- glmFit(dge_sub, design)

      contrast <- makeContrasts(
        contrasts = paste0("`", cond1, "` - `", cond2, "`"),
        levels = design
      )
      lrt <- glmLRT(fit, contrast = contrast)

      res <- topTags(lrt, n = Inf)$table
      res$Trans <- rownames(res)
      sig <- res[which(res$FDR < 0.05), ]

      if (nrow(sig) > 0) {
        UP <- subset(sig, logFC >= 0)
        DOWN <- subset(sig, logFC < 0)

        # Format
        UP_formatted <- UP %>%
          select(Trans, logFC, logCPM, LR, PValue, FDR) %>%
          mutate(logCPM = log2(logCPM + 1))
        DOWN_formatted <- DOWN %>%
          select(Trans, logFC, logCPM, LR, PValue, FDR) %>%
          mutate(logCPM = log2(logCPM + 1))

        # Annotate
        UP_anno <- left_join(UP_formatted, NorCount_anno, by = "Trans")
        DOWN_anno <- left_join(DOWN_formatted, NorCount_anno, by = "Trans")

        # Save to Excel
        up_file <- file.path(output_dir, "UP_edgeR.xlsx")
        down_file <- file.path(output_dir, "DOWN_edgeR.xlsx")

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
write.xlsx(results_summary, file.path(output_dir, "DEG_Analysis_Summary_edgeR.xlsx"), rowNames = FALSE)
message("=== EDGE-R ANALYSIS COMPLETE ===")
print(results_summary)
