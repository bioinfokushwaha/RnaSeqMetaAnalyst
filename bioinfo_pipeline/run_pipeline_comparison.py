#!/usr/bin/env python3
"""
RNA-Seq Pipeline Correlation Analysis
=====================================
This script performs comprehensive analysis of pipeline combinations:
1. Loads normalized count matrices from 16 combinations
2. Identifies gene ID type and merges with core_genes.xlsx
3. Counts non-zero genes across samples
4. Generates cross-comparison tables for UP/DOWN regulated genes
5. Creates lower triangular correlation matrices with heatmaps
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Pipeline combination codes (3-letter)
# Format: Aligner(1) + Quantifier(1) + DEG_tool(1)
PIPELINE_CODES = {
    'BFE': {'aligner': 'Bowtie2', 'quant': 'FC', 'deg': 'edgeR'},
    'BFD': {'aligner': 'Bowtie2', 'quant': 'FC', 'deg': 'DESeq2'},
    'BHE': {'aligner': 'Bowtie2', 'quant': 'HTSeq', 'deg': 'edgeR'},
    'BHD': {'aligner': 'Bowtie2', 'quant': 'HTSeq', 'deg': 'DESeq2'},
    'BRE': {'aligner': 'Bowtie2', 'quant': 'RSEM', 'deg': 'edgeR'},
    'BRD': {'aligner': 'Bowtie2', 'quant': 'RSEM', 'deg': 'DESeq2'},
    'HFE': {'aligner': 'HISAT2', 'quant': 'FC', 'deg': 'edgeR'},
    'HFD': {'aligner': 'HISAT2', 'quant': 'FC', 'deg': 'DESeq2'},
    'HHE': {'aligner': 'HISAT2', 'quant': 'HTSeq', 'deg': 'edgeR'},
    'HHD': {'aligner': 'HISAT2', 'quant': 'HTSeq', 'deg': 'DESeq2'},
    'SFE': {'aligner': 'STAR', 'quant': 'FC', 'deg': 'edgeR'},
    'SFD': {'aligner': 'STAR', 'quant': 'FC', 'deg': 'DESeq2'},
    'SHE': {'aligner': 'STAR', 'quant': 'HTSeq', 'deg': 'edgeR'},
    'SHD': {'aligner': 'STAR', 'quant': 'HTSeq', 'deg': 'DESeq2'},
    'SRE': {'aligner': 'STAR', 'quant': 'RSEM', 'deg': 'edgeR'},
    'SRD': {'aligner': 'STAR', 'quant': 'RSEM', 'deg': 'DESeq2'},
}

# File patterns for normalized count matrices
FILE_PATTERNS = {
    'BFE': 'B_FC_counts_clean_edgeR_normalized.csv',
    'BFD': 'B_FC_counts_clean_DESeq2_normalized.csv',
    'BHE': 'B_HTSeq_counts_clean_edgeR_normalized.csv',
    'BHD': 'B_HTSeq_counts_clean_DESeq2_normalized.csv',
    'BRE': 'B_RSEM_counts_clean_edgeR_normalized.csv',
    'BRD': 'B_RSEM_counts_clean_DESeq2_normalized.csv',
    'HFE': 'H_FC_counts_clean_edgeR_normalized.csv',
    'HFD': 'H_FC_counts_clean_DESeq2_normalized.csv',
    'HHE': 'H_HTSeq_counts_clean_edgeR_normalized.csv',
    'HHD': 'H_HTSeq_counts_clean_DESeq2_normalized.csv',
    'SFE': 'S_FC_counts_clean_edgeR_normalized.csv',
    'SFD': 'S_FC_counts_clean_DESeq2_normalized.csv',
    'SHE': 'S_HTSeq_counts_clean_edgeR_normalized.csv',
    'SHD': 'S_HTSeq_counts_clean_DESeq2_normalized.csv',
    'SRE': 'S_RSEM_counts_clean_edgeR_normalized.csv',
    'SRD': 'S_RSEM_counts_clean_DESeq2_normalized.csv',
}

# DEG result file patterns
DEG_FILE_PATTERNS = {
    'UP': '{}_UP.xlsx',
    'DOWN': '{}_DOWN.xlsx',
}


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def detect_gene_id_type(gene_id):
    """
    Detect the type of gene identifier
    
    Returns: 'ensembl_gene', 'gene_name', 'gene_id', or 'transcript_id'
    """
    gene_id = str(gene_id).strip()
    
    if gene_id.startswith('ENS'):
        if 'T' in gene_id and gene_id.index('T') < 15:  # ENST = transcript
            return 'transcript_id'
        else:  # ENSG = gene
            return 'ensembl_gene'
    elif gene_id.startswith('gene:'):
        return 'gene_id'
    elif gene_id.startswith('LOC') or gene_id.startswith('NM_') or gene_id.startswith('XM_'):
        return 'transcript_id'
    else:
        # Assume gene name if alphabetic
        if gene_id.replace('-', '').replace('_', '').isalpha():
            return 'gene_name'
        else:
            return 'gene_id'


def detect_id_column(df, core_genes):
    """
    Automatically detect which column type to use for merging
    
    Returns: (column_name_in_normalized, column_name_in_core_genes)
    """
    # Get first gene ID from normalized file (skip header)
    first_id = None
    for col in df.columns:
        if col.lower() in ['geneid', 'gene', 'gene_id', 'id']:
            first_id = str(df[col].iloc[0])
            break
    
    if first_id is None:
        first_id = str(df.iloc[0, 0])  # Use first column
    
    # Detect type
    id_type = detect_gene_id_type(first_id)
    
    # Map to core_genes column
    type_mapping = {
        'ensembl_gene': 'ensembl_geneid',
        'gene_name': 'gene_name',
        'gene_id': 'gene_id',
        'transcript_id': 'transcript_id'
    }
    
    # Find normalized file column
    norm_col = None
    for col in df.columns:
        if col.lower() in ['geneid', 'gene', 'gene_id', 'id', 'gene_name']:
            norm_col = col
            break
    
    if norm_col is None:
        norm_col = df.columns[0]
    
    # Find core_genes column
    core_col = None
    target_col = type_mapping.get(id_type, 'gene_id')
    
    for col in core_genes.columns:
        if col.lower().replace('_', '').replace(' ', '') == target_col.lower().replace('_', '').replace(' ', ''):
            core_col = col
            break
    
    if core_col is None:
        # Fallback: try to match any column
        for col in core_genes.columns:
            core_col = col
            break
    
    print(f"  Detected ID type: {id_type}")
    print(f"  Normalized file column: {norm_col}")
    print(f"  Core genes column: {core_col}")
    
    return norm_col, core_col


def load_normalized_counts(deg_dir, pipeline_code):
    """Load normalized count matrix for a pipeline combination"""
    file_pattern = FILE_PATTERNS.get(pipeline_code)
    if not file_pattern:
        return None
    
    file_path = deg_dir / file_pattern
    
    if not file_path.exists():
        print(f"  ⚠ File not found: {file_path}")
        return None
    
    try:
        df = pd.read_csv(file_path)
        print(f"  ✓ Loaded {pipeline_code}: {df.shape}")
        return df
    except Exception as e:
        print(f"  ✗ Error loading {pipeline_code}: {e}")
        return None


def load_deg_results(deg_dir, pipeline_code, direction):
    """Load DEG results (UP or DOWN regulated genes)"""
    # Construct filename based on pipeline code
    info = PIPELINE_CODES[pipeline_code]
    aligner_prefix = pipeline_code[0]  # B, H, or S
    
    # Map quantifier
    quant_map = {'FC': 'FC', 'HTSeq': 'HTSeq', 'RSEM': 'RSEM'}
    quant = quant_map.get(info['quant'], info['quant'])
    
    # Construct filename pattern: B_FC_counts_clean_edgeR_UP.xlsx
    filename = f"{aligner_prefix}_{quant}_counts_clean_{info['deg']}_{direction}.xlsx"
    file_path = deg_dir / filename
    
    if not file_path.exists():
        print(f"  ⚠ DEG file not found: {file_path}")
        return None
    
    try:
        df = pd.read_excel(file_path)
        print(f"  ✓ Loaded {pipeline_code} {direction}: {len(df)} genes")
        return df
    except Exception as e:
        print(f"  ✗ Error loading {pipeline_code} {direction}: {e}")
        return None


def count_nonzero_genes(df, gene_col):
    """Count non-zero genes for each sample"""
    # Get sample columns (exclude gene ID column)
    sample_cols = [col for col in df.columns if col != gene_col]
    
    nonzero_counts = {}
    for col in sample_cols:
        try:
            nonzero = (df[col] > 0).sum()
            nonzero_counts[col] = nonzero
        except:
            nonzero_counts[col] = 0
    
    return nonzero_counts


def create_comparison_matrix(deg_dir, direction='UP'):
    """Create 16x16 comparison matrix for UP or DOWN regulated genes"""
    print(f"\n>>> Creating {direction}-regulated gene comparison matrix...")
    
    # Initialize matrix
    codes = list(PIPELINE_CODES.keys())
    matrix = pd.DataFrame(index=codes, columns=codes, dtype=int)
    
    # Load all DEG results
    deg_data = {}
    for code in codes:
        deg_df = load_deg_results(deg_dir, code, direction)
        if deg_df is not None:
            # Get gene IDs (first column usually)
            gene_col = deg_df.columns[0]
            deg_data[code] = set(deg_df[gene_col].astype(str))
    
    # Calculate overlaps
    for i, code1 in enumerate(codes):
        for j, code2 in enumerate(codes):
            if code1 in deg_data and code2 in deg_data:
                overlap = len(deg_data[code1] & deg_data[code2])
                matrix.loc[code1, code2] = overlap
            else:
                matrix.loc[code1, code2] = 0
    
    return matrix.astype(int)


def plot_lower_triangular_heatmap(matrix, title, output_path):
    """Plot lower triangular heatmap with correlation values"""
    print(f"\n>>> Creating heatmap: {title}")
    
    # Convert to numpy array
    data = matrix.values.astype(float)
    
    # Create mask for upper triangle
    mask = np.triu(np.ones_like(data, dtype=bool))
    
    # Create figure
    plt.figure(figsize=(14, 12))
    
    # Plot heatmap
    sns.heatmap(
        data,
        mask=mask,
        annot=True,
        fmt='g',
        cmap='coolwarm',
        center=data.mean(),
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8, "label": "Number of Common Genes"},
        xticklabels=matrix.columns,
        yticklabels=matrix.index
    )
    
    plt.title(title, fontsize=16, pad=20)
    plt.xlabel('Pipeline Combination', fontsize=12)
    plt.ylabel('Pipeline Combination', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_path}")
    plt.close()


def plot_correlation_heatmap(matrix, title, output_path):
    """Plot correlation heatmap (lower triangular)"""
    print(f"\n>>> Creating correlation heatmap: {title}")
    
    # Calculate correlation
    corr_matrix = matrix.corr()
    
    # Create mask for upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    
    # Create figure
    plt.figure(figsize=(14, 12))
    
    # Plot heatmap
    sns.heatmap(
        corr_matrix,
        mask=mask,
        annot=True,
        fmt='.2f',
        cmap='coolwarm',
        vmax=1,
        vmin=-1,
        center=0,
        square=True,
        linewidths=0.5,
        cbar_kws={"shrink": 0.8, "label": "Correlation Coefficient"},
        xticklabels=corr_matrix.columns,
        yticklabels=corr_matrix.index
    )
    
    plt.title(title, fontsize=16, pad=20)
    plt.xlabel('Pipeline Combination', fontsize=12)
    plt.ylabel('Pipeline Combination', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_path}")
    plt.close()


# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================

def main():
    """Main analysis pipeline"""
    
    print("")
    print("=" * 70)
    print("RNA-SEQ PIPELINE CORRELATION ANALYSIS")
    print("=" * 70)
    print("")
    
    # Setup paths
    deg_dir = Path('/data/results/DEG')
    output_dir = Path('/data/results/pipeline_comparison')
    output_dir.mkdir(exist_ok=True)
    
    core_genes_path = Path('/data/data/core_genes.xlsx')
    
    # Check if core_genes.xlsx exists
    if not core_genes_path.exists():
        print(f"⚠ WARNING: core_genes.xlsx not found at {core_genes_path}")
        print("  Proceeding without gene annotation...")
        core_genes = None
    else:
        print(f"✓ Loading core genes from: {core_genes_path}")
        core_genes = pd.read_excel(core_genes_path)
        print(f"  Loaded {len(core_genes)} genes with columns: {list(core_genes.columns)}")
    
    # ========================================================================
    # PART 1: NON-ZERO GENE COUNTS ANALYSIS
    # ========================================================================
    print("\n" + "=" * 70)
    print("PART 1: NON-ZERO GENE COUNTS ANALYSIS")
    print("=" * 70)
    
    nonzero_summary = {}
    merged_data = {}
    
    for code in PIPELINE_CODES.keys():
        print(f"\n>>> Processing {code}: {PIPELINE_CODES[code]}")
        
        # Load normalized counts
        norm_df = load_normalized_counts(deg_dir, code)
        if norm_df is None:
            continue
        
        # Detect gene ID column
        gene_col = None
        for col in norm_df.columns:
            if col.lower() in ['geneid', 'gene', 'gene_id', 'id']:
                gene_col = col
                break
        if gene_col is None:
            gene_col = norm_df.columns[0]
        
        # Merge with core_genes if available
        if core_genes is not None:
            try:
                norm_col, core_col = detect_id_column(norm_df, core_genes)
                merged = pd.merge(
                    norm_df,
                    core_genes,
                    left_on=norm_col,
                    right_on=core_col,
                    how='inner'
                )
                print(f"  ✓ Merged: {len(merged)} genes matched")
                merged_data[code] = merged
                
                # Count non-zeros on merged data
                nonzero_counts = count_nonzero_genes(merged, norm_col)
            except Exception as e:
                print(f"  ⚠ Merge failed: {e}")
                nonzero_counts = count_nonzero_genes(norm_df, gene_col)
        else:
            nonzero_counts = count_nonzero_genes(norm_df, gene_col)
        
        nonzero_summary[code] = nonzero_counts
        print(f"  ✓ Non-zero counts calculated for {len(nonzero_counts)} samples")
    
    # Save non-zero counts summary
    if nonzero_summary:
        print("\n>>> Saving non-zero gene counts summary...")
        nonzero_df = pd.DataFrame(nonzero_summary).T
        nonzero_df.index.name = 'Pipeline'
        output_file = output_dir / 'nonzero_gene_counts_by_sample.xlsx'
        nonzero_df.to_excel(output_file)
        print(f"  ✓ Saved: {output_file}")
        
        # Also save as CSV for easy viewing
        csv_file = output_dir / 'nonzero_gene_counts_by_sample.csv'
        nonzero_df.to_csv(csv_file)
        print(f"  ✓ Saved: {csv_file}")
    
    # ========================================================================
    # PART 2: UP-REGULATED GENES COMPARISON
    # ========================================================================
    print("\n" + "=" * 70)
    print("PART 2: UP-REGULATED GENES COMPARISON")
    print("=" * 70)
    
    up_matrix = create_comparison_matrix(deg_dir, 'UP')
    
    # Save matrix
    output_file = output_dir / 'upregulated_genes_comparison_matrix.xlsx'
    up_matrix.to_excel(output_file)
    print(f"\n✓ Saved UP-regulated comparison matrix: {output_file}")
    
    # Plot heatmap
    plot_lower_triangular_heatmap(
        up_matrix,
        'UP-Regulated Genes: Pipeline Comparison (Lower Triangle)',
        output_dir / 'upregulated_genes_heatmap.png'
    )
    
    # ========================================================================
    # PART 3: DOWN-REGULATED GENES COMPARISON
    # ========================================================================
    print("\n" + "=" * 70)
    print("PART 3: DOWN-REGULATED GENES COMPARISON")
    print("=" * 70)
    
    down_matrix = create_comparison_matrix(deg_dir, 'DOWN')
    
    # Save matrix
    output_file = output_dir / 'downregulated_genes_comparison_matrix.xlsx'
    down_matrix.to_excel(output_file)
    print(f"\n✓ Saved DOWN-regulated comparison matrix: {output_file}")
    
    # Plot heatmap
    plot_lower_triangular_heatmap(
        down_matrix,
        'DOWN-Regulated Genes: Pipeline Comparison (Lower Triangle)',
        output_dir / 'downregulated_genes_heatmap.png'
    )
    
    # ========================================================================
    # PART 4: CORRELATION ANALYSIS
    # ========================================================================
    print("\n" + "=" * 70)
    print("PART 4: CORRELATION ANALYSIS")
    print("=" * 70)
    
    # Plot correlation heatmap for UP-regulated genes
    plot_correlation_heatmap(
        up_matrix,
        'UP-Regulated Genes: Correlation Matrix (Lower Triangle)',
        output_dir / 'upregulated_genes_correlation.png'
    )
    
    # Plot correlation heatmap for DOWN-regulated genes
    plot_correlation_heatmap(
        down_matrix,
        'DOWN-Regulated Genes: Correlation Matrix (Lower Triangle)',
        output_dir / 'downregulated_genes_correlation.png'
    )
    
    # ========================================================================
    # PART 5: SUMMARY REPORT
    # ========================================================================
    print("\n" + "=" * 70)
    print("GENERATING SUMMARY REPORT")
    print("=" * 70)
    
    summary_lines = [
        "RNA-Seq Pipeline Comparison Analysis Summary",
        "=" * 60,
        "",
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"DEG Directory: {deg_dir}",
        f"Output Directory: {output_dir}",
        "",
        "Pipeline Combinations Analyzed:",
        "-" * 60,
    ]
    
    for code, info in PIPELINE_CODES.items():
        summary_lines.append(f"  {code}: {info['aligner']}-{info['quant']}-{info['deg']}")
    
    summary_lines.extend([
        "",
        "Results Generated:",
        "-" * 60,
        f"  ✓ Non-zero gene counts by sample",
        f"  ✓ UP-regulated genes comparison (16x16 matrix)",
        f"  ✓ DOWN-regulated genes comparison (16x16 matrix)",
        f"  ✓ Correlation matrices (lower triangular)",
        f"  ✓ Heatmaps for all comparisons",
        "",
        "Files Created:",
        "-" * 60,
        "  • nonzero_gene_counts_by_sample.xlsx",
        "  • nonzero_gene_counts_by_sample.csv",
        "  • upregulated_genes_comparison_matrix.xlsx",
        "  • upregulated_genes_heatmap.png",
        "  • upregulated_genes_correlation.png",
        "  • downregulated_genes_comparison_matrix.xlsx",
        "  • downregulated_genes_heatmap.png",
        "  • downregulated_genes_correlation.png",
        "",
    ])
    
    # Add statistics
    if up_matrix is not None:
        summary_lines.extend([
            "UP-Regulated Genes Statistics:",
            "-" * 60,
            f"  Mean overlap: {up_matrix.values[~np.eye(16, dtype=bool)].mean():.0f} genes",
            f"  Max overlap: {up_matrix.values[~np.eye(16, dtype=bool)].max():.0f} genes",
            f"  Min overlap: {up_matrix.values[~np.eye(16, dtype=bool)].min():.0f} genes",
            "",
        ])
    
    if down_matrix is not None:
        summary_lines.extend([
            "DOWN-Regulated Genes Statistics:",
            "-" * 60,
            f"  Mean overlap: {down_matrix.values[~np.eye(16, dtype=bool)].mean():.0f} genes",
            f"  Max overlap: {down_matrix.values[~np.eye(16, dtype=bool)].max():.0f} genes",
            f"  Min overlap: {down_matrix.values[~np.eye(16, dtype=bool)].min():.0f} genes",
            "",
        ])
    
    summary_lines.extend([
        "=" * 60,
        "Analysis Complete!",
        ""
    ])
    
    summary_text = "\n".join(summary_lines)
    
    # Print summary
    print("\n" + summary_text)
    
    # Save summary
    summary_file = output_dir / 'analysis_summary.txt'
    with open(summary_file, 'w') as f:
        f.write(summary_text)
    print(f"✓ Summary saved: {summary_file}")
    
    print("\n" + "=" * 70)
    print("✅ PIPELINE COMPARISON ANALYSIS COMPLETE!")
    print("=" * 70)
    print(f"\n📁 Results location: {output_dir}")
    print("")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠ Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
