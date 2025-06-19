# Load libraries
library(edgeR)
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

# Create DGEList object
dge <- DGEList(counts = count_matrix, group = sampleinfo$CONDITION)

# Normalize
dge <- calcNormFactors(dge)

# Normalized counts for annotation
cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE)
NorCount_anno <- data.frame(
  Trans = rownames(cpm_norm),
  NormalizedCount_Mean = rowMeans(cpm_norm),
  stringsAsFactors = FALSE
)

# Create output directory
dir.create("DEG_Results_edgeR", showWarnings = FALSE)

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
        # Subset data for current pair
        keep_samples <- sampleinfo$CONDITION %in% c(cond1, cond2)
        dge_sub <- dge[, keep_samples]
        dge_sub$samples$group <- droplevels(factor(sampleinfo$CONDITION[keep_samples]))
        
        # Syntactically valid names
        group_clean <- make.names(dge_sub$samples$group)
        design <- model.matrix(~0 + group_clean)
        colnames(design) <- levels(factor(group_clean))
        
        # Fit model
        dge_sub <- estimateDisp(dge_sub, design)
        fit <- glmFit(dge_sub, design)
        
        # Valid contrast names
        cond1_valid <- make.names(cond1)
        cond2_valid <- make.names(cond2)
        contrast <- makeContrasts(contrasts = paste0(cond1_valid, "-", cond2_valid), levels = design)
        lrt <- glmLRT(fit, contrast = contrast)
        
        # Filter DEGs
        res <- topTags(lrt, n = Inf)$table
        sig <- res[which(res$FDR < 0.05), ]
        print(paste("Total significant DEGs:", nrow(sig)))
        
        if (nrow(sig) > 0) {
          # UP and DOWN
          UP <- subset(sig, logFC >= 0)
          UP$Trans <- rownames(UP)
          DOWN <- subset(sig, logFC < 0)
          DOWN$Trans <- rownames(DOWN)
          
          print(paste("UP:", nrow(UP), "DOWN:", nrow(DOWN)))
          
          # Annotate
          UP_anno <- left_join(UP, NorCount_anno, by = "Trans")
          DOWN_anno <- left_join(DOWN, NorCount_anno, by = "Trans")
          
          # Save UP
          up_file <- "DEG_Results_edgeR/Treatment_Rgene_UP_edgeR.xlsx"
          wb_up <- if (file.exists(up_file)) loadWorkbook(up_file) else createWorkbook()
          addWorksheet(wb_up, comp_name)
          writeData(wb_up, comp_name, UP_anno, rowNames = FALSE)
          saveWorkbook(wb_up, up_file, overwrite = TRUE)
          
          # Save DOWN
          down_file <- "DEG_Results_edgeR/Treatment_Rgene_DOWN_edgeR.xlsx"
          wb_down <- if (file.exists(down_file)) loadWorkbook(down_file) else createWorkbook()
          addWorksheet(wb_down, comp_name)
          writeData(wb_down, comp_name, DOWN_anno, rowNames = FALSE)
          saveWorkbook(wb_down, down_file, overwrite = TRUE)
          
          # Track
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
write.xlsx(results_summary, "DEG_Results_edgeR/DEG_Analysis_Summary_edgeR.xlsx", rowNames = FALSE)

# Print summary
print("=== EDGE-R ANALYSIS COMPLETE ===")
print(results_summary)
