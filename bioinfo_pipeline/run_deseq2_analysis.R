#!/usr/bin/env Rscript

# ===========================
# DESeq2 Differential Expression Analysis
# ===========================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
})

# ===========================
set.seed(0)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript run_DESeq2_analysis.R <count_file> <sampleinfo_file>")
}

count_file <- args[1]
sampleinfo_file <- args[2]

LOG2FC_THRESHOLD <- 1
FDR_THRESHOLD <- 0.05

# ===========================
filename <- tools::file_path_sans_ext(basename(count_file))
output_dir <- file.path("results/DEG", filename)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================
counts <- read.csv(count_file, row.names=1, check.names=FALSE)
sampleinfo <- read.delim(sampleinfo_file, sep="\t")
colnames(sampleinfo) <- trimws(colnames(sampleinfo))
sampleinfo$Sample <- trimws(sampleinfo$Sample)
sampleinfo$Condition <- trimws(sampleinfo$Condition)

# ===========================
# Match samples
count_cols_clean <- gsub("\\.bam$|_Aligned.*$|\\.fastq.*$|\\.htseq.*$", "", colnames(counts))
sampleinfo$CleanName <- gsub("\\.bam$|_Aligned.*$|\\.fastq.*$|\\.htseq.*$", "", sampleinfo$Sample)
matching <- intersect(count_cols_clean, sampleinfo$CleanName)
if(length(matching)==0) stop("No matching samples found!")
count_idx <- match(sampleinfo$CleanName, count_cols_clean)
valid_idx <- !is.na(count_idx)
count_matrix <- round(as.matrix(counts[, count_idx[valid_idx]]))
sampleinfo_matched <- sampleinfo[valid_idx, ]

# ===========================
# Filter housekeeping / non-zero
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

# ===========================
colData <- data.frame(CONDITION=factor(sampleinfo_matched$Condition), row.names=colnames(count_matrix))
dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=colData, design=~CONDITION)

dds <- DESeq(dds)
norm_counts <- counts(dds, normalized=TRUE)
write.xlsx(data.frame(Gene=rownames(norm_counts), norm_counts), file.path(output_dir,paste0(filename,"_DESeq2_Normalized.xlsx")), rowNames=FALSE)

# ===========================
# PCA plot
pca <- prcomp(t(norm_counts))
var_pct <- round(100*(pca$sdev^2)/sum(pca$sdev^2),2)
scores <- data.frame(pca$x, Condition=colData$CONDITION)
write.xlsx(scores, file.path(output_dir,paste0(filename,"_DESeq2_PCA_Scores.xlsx")), rowNames=FALSE)

p <- ggplot(scores, aes(PC1, PC2, color=Condition, label=rownames(scores))) +
  geom_point(size=4) + geom_text(vjust=-0.8, hjust=0.5, size=3, show.legend=FALSE) +
  labs(x=paste0("PC1 (",var_pct[1],"%)"), y=paste0("PC2 (",var_pct[2],"%)"),
       title=paste("DESeq2 PCA -",filename)) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=14))
ggsave(file.path(output_dir,paste0(filename,"_DESeq2_PCA_Plot.png")), p, width=10, height=8, dpi=300)

# ===========================
# Differential Expression per pairwise comparison
results_summary <- data.frame(Comparison=character(), UP=numeric(), DOWN=numeric(), Total_DEGs=numeric())
conditions <- levels(colData$CONDITION)

for(i in 1:(length(conditions)-1)){
  for(j in (i+1):length(conditions)){
    c1 <- conditions[i]; c2 <- conditions[j]; comp <- paste0(c1,"_vs_",c2)
    res <- results(dds, contrast=c("CONDITION", c1, c2))
    res_df <- as.data.frame(res)
    res_df$Gene <- rownames(res_df)
    
    # Volcano info
    res_df$Status <- "NS"
    res_df$Status[res_df$padj<FDR_THRESHOLD & res_df$log2FoldChange>LOG2FC_THRESHOLD] <- "UP"
    res_df$Status[res_df$padj<FDR_THRESHOLD & res_df$log2FoldChange< -LOG2FC_THRESHOLD] <- "DOWN"
    
    # Volcano plot
    png(file.path(output_dir,paste0(filename,"_",comp,"_DESeq2_Volcano.png")), width=3400, height=3400, res=600)
    print(
      ggplot(res_df, aes(log2FoldChange, -log10(padj), color=Status)) +
        geom_point(alpha=0.6, size=2) +
        geom_vline(xintercept=c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype="dashed") +
        geom_hline(yintercept=-log10(FDR_THRESHOLD), linetype="dashed") +
        scale_color_manual(values=c(UP="red", DOWN="blue", NS="grey")) +
        theme_bw() +
        labs(title=paste("Volcano:",filename,"-",comp))
    )
    dev.off()
    
    # Save UP/DOWN
    UP <- res_df[res_df$Status=="UP",]; DOWN <- res_df[res_df$Status=="DOWN",]
    if(nrow(UP)>0) write.xlsx(UP, file.path(output_dir,paste0(filename,"_",comp,"_DESeq2_UP.xlsx")), rowNames=FALSE)
    if(nrow(DOWN)>0) write.xlsx(DOWN, file.path(output_dir,paste0(filename,"_",comp,"_DESeq2_DOWN.xlsx")), rowNames=FALSE)
    
    results_summary <- rbind(results_summary, data.frame(Comparison=comp, UP=nrow(UP), DOWN=nrow(DOWN), Total_DEGs=nrow(res_df[res_df$padj<FDR_THRESHOLD,])))
  }
}

write.xlsx(results_summary, file.path(output_dir,paste0(filename,"_DESeq2_DEG_Summary.xlsx")), rowNames=FALSE)
cat("DESeq2 Analysis complete! Output in:", output_dir,"\n")
