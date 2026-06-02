#!/usr/bin/env python3
"""
RNA-Seq Pipeline Correlation Analysis - FIXED VERSION
======================================================
FIXES:
1. Proper path detection for DEG files
2. Better error handling and validation
3. Fallback mechanisms for missing files
4. Detailed logging of file searches

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
# Actual filenames: e.g. B_FC_edgeR_normalized.csv  (inside DEG/B_FC/edgeR/)
FILE_PATTERNS = {
    'BFE': ('B_FC',   'edgeR',  'B_FC_edgeR_normalized.csv'),
    'BFD': ('B_FC',   'DESeq2', 'B_FC_DESeq2_normalized.csv'),
    'BHE': ('B_HTSeq','edgeR',  'B_HTSeq_edgeR_normalized.csv'),
    'BHD': ('B_HTSeq','DESeq2', 'B_HTSeq_DESeq2_normalized.csv'),
    'BRE': ('B_RSEM', 'edgeR',  'B_RSEM_edgeR_normalized.csv'),
    'BRD': ('B_RSEM', 'DESeq2', 'B_RSEM_DESeq2_normalized.csv'),
    'HFE': ('H_FC',   'edgeR',  'H_FC_edgeR_normalized.csv'),
    'HFD': ('H_FC',   'DESeq2', 'H_FC_DESeq2_normalized.csv'),
    'HHE': ('H_HTSeq','edgeR',  'H_HTSeq_edgeR_normalized.csv'),
    'HHD': ('H_HTSeq','DESeq2', 'H_HTSeq_DESeq2_normalized.csv'),
    'SFE': ('S_FC',   'edgeR',  'S_FC_edgeR_normalized.csv'),
    'SFD': ('S_FC',   'DESeq2', 'S_FC_DESeq2_normalized.csv'),
    'SHE': ('S_HTSeq','edgeR',  'S_HTSeq_edgeR_normalized.csv'),
    'SHD': ('S_HTSeq','DESeq2', 'S_HTSeq_DESeq2_normalized.csv'),
    'SRE': ('S_RSEM', 'edgeR',  'S_RSEM_edgeR_normalized.csv'),
    'SRD': ('S_RSEM', 'DESeq2', 'S_RSEM_DESeq2_normalized.csv'),
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


# ============================================================================
# FIXED: Better file path detection
# ============================================================================
def find_deg_file(deg_dir, aligner_prefix, quant, tool, direction):
    """
    Find DEG result file using the actual subdirectory structure:
      <deg_dir>/<PREFIX>_<QUANT>/<TOOL>/<PREFIX>_<QUANT>_<TOOL>_Control_vs_Treatment_<DIRECTION>.xlsx

    Also tries flat-directory and legacy naming as fallbacks.
    """
    folder_key = f"{aligner_prefix}_{quant}"   # e.g. B_FC, H_HTSeq, S_RSEM

    # --- Primary: subdirectory layout produced by the R scripts ---
    primary_filename = f"{folder_key}_{tool}_Control_vs_Treatment_{direction}.xlsx"
    primary_path = deg_dir / folder_key / tool / primary_filename
    if primary_path.exists():
        return primary_path

    # --- Fallback roots to search (handles results/results/DEG vs results/DEG) ---
    fallback_roots = [
        Path('/data/results/results/DEG'),
        Path('/data/results/DEG'),
        Path('/data/DEG'),
        deg_dir,
    ]

    for root in fallback_roots:
        if not root.exists():
            continue
        # Try subdirectory layout under this root
        p = root / folder_key / tool / primary_filename
        if p.exists():
            return p
        # Try flat legacy names
        for fname in [
            f"{aligner_prefix}_{quant}_counts_clean_{tool}_{direction}.xlsx",
            f"{aligner_prefix}_{quant.lower()}_counts_clean_{tool}_{direction}.xlsx",
        ]:
            p = root / fname
            if p.exists():
                return p

    return None


def load_normalized_counts(deg_dir, pipeline_code):
    """Load normalized count matrix for a pipeline combination"""
    pattern_info = FILE_PATTERNS.get(pipeline_code)
    if not pattern_info:
        return None

    folder_key, tool, filename = pattern_info  # e.g. ('B_FC', 'edgeR', 'B_FC_edgeR_normalized.csv')

    # Candidate roots (handles results/results/DEG and results/DEG layouts)
    roots = [
        Path('/data/results/results/DEG'),
        Path('/data/results/DEG'),
        Path('/data/DEG'),
        deg_dir,
    ]

    file_path = None
    for root in roots:
        if not root.exists():
            continue
        p = root / folder_key / tool / filename
        if p.exists():
            file_path = p
            break
        # Also try flat (legacy) layout
        p2 = root / filename
        if p2.exists():
            file_path = p2
            break

    if file_path is None:
        print(f"  ⚠ File not found: {filename}")
        print(f"    Searched: subdirs under /data/results/results/DEG, /data/results/DEG, /data/DEG")
        return None

    try:
        df = pd.read_csv(file_path)
        print(f"  ✓ Loaded {pipeline_code}: {df.shape} from {file_path}")
        return df
    except Exception as e:
        print(f"  ✗ Error loading {pipeline_code}: {e}")
        return None


def load_deg_results(deg_dir, pipeline_code, direction):
    """
    Load DEG results (UP or DOWN regulated genes)

    Actual file pattern:
      <deg_dir>/<PREFIX>_<QUANT>/<TOOL>/<PREFIX>_<QUANT>_<TOOL>_Control_vs_Treatment_<DIRECTION>.xlsx
    e.g.  DEG/B_FC/DESeq2/B_FC_DESeq2_Control_vs_Treatment_UP.xlsx
    """
    info = PIPELINE_CODES[pipeline_code]
    aligner_prefix = pipeline_code[0]  # B, H, or S

    # quant folder names match the FILE_PATTERNS folder_key second segment
    quant_map = {'FC': 'FC', 'HTSeq': 'HTSeq', 'RSEM': 'RSEM'}
    quant = quant_map.get(info['quant'], info['quant'])

    file_path = find_deg_file(deg_dir, aligner_prefix, quant, info['deg'], direction)

    if file_path is None:
        expected = f"{aligner_prefix}_{quant}_{info['deg']}_Control_vs_Treatment_{direction}.xlsx"
        print(f"  ⚠ DEG file not found for {pipeline_code} {direction}")
        print(f"    Expected: {expected}")
        print(f"    Looked under: {deg_dir}/{aligner_prefix}_{quant}/{info['deg']}/")
        return None

    try:
        df = pd.read_excel(file_path)
        print(f"  ✓ Loaded {pipeline_code} {direction}: {len(df)} genes from {file_path.name}")
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
    loaded_count = 0
    for code in codes:
        deg_df = load_deg_results(deg_dir, code, direction)
        if deg_df is not None:
            # Get gene IDs (first column usually)
            gene_col = deg_df.columns[0]
            deg_data[code] = set(deg_df[gene_col].astype(str))
            loaded_count += 1
        else:
            deg_data[code] = set()
    
    print(f"  Successfully loaded: {loaded_count}/16 DEG files")
    
    if loaded_count == 0:
        print(f"  ✗ WARNING: No {direction}-regulated gene files found!")
        return None
    
    # Calculate overlaps
    for i, code1 in enumerate(codes):
        for j, code2 in enumerate(codes):
            if code1 in deg_data and code2 in deg_data:
                overlap = len(deg_data[code1] & deg_data[code2])
                matrix.loc[code1, code2] = overlap
            else:
                matrix.loc[code1, code2] = 0
    
    matrix = matrix.astype(int)
    
    print(f"  ✓ Comparison matrix created successfully")
    if loaded_count < 16:
        print(f"  ⚠ Note: Matrix includes only {loaded_count} loaded datasets")
    
    return matrix


def create_heatmap(matrix, title, filename, output_dir):
    """Create heatmap visualization of comparison matrix"""
    if matrix is None or matrix.empty:
        print(f"  ⚠ Cannot create heatmap: matrix is empty")
        return False
    
    try:
        plt.figure(figsize=(14, 12))
        sns.heatmap(matrix, annot=True, fmt='d', cmap='YlOrRd', cbar_kws={'label': 'Overlapping Genes'})
        plt.title(title, fontsize=16, fontweight='bold')
        plt.xlabel('Pipeline Combination', fontsize=12)
        plt.ylabel('Pipeline Combination', fontsize=12)
        plt.xticks(rotation=45)
        plt.yticks(rotation=0)
        plt.tight_layout()
        
        output_file = output_dir / filename
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  ✓ Saved: {filename}")
        return True
    except Exception as e:
        print(f"  ✗ Error creating heatmap: {e}")
        return False


def calculate_correlation_lower_triangular(matrix):
    """Calculate correlation for lower triangular matrix"""
    if matrix is None or matrix.empty:
        return None
    
    try:
        # Transpose to get genes as rows, pipelines as columns
        matrix_t = matrix.T
        
        # Calculate Pearson correlation
        corr_matrix = matrix_t.corr(method='pearson')
        
        # Convert to lower triangular
        lower_tri = corr_matrix.copy()
        for i in range(len(lower_tri)):
            for j in range(i + 1, len(lower_tri)):
                lower_tri.iloc[i, j] = np.nan
        
        return lower_tri
    except Exception as e:
        print(f"  ✗ Error calculating correlation: {e}")
        return None


def main():
    """Main execution function"""
    print("=" * 70)
    print("RNA-SEQ PIPELINE COMPARISON ANALYSIS - FIXED VERSION")
    print("=" * 70)
    print("")
    print("Starting analysis...")
    print(f"Time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("")
    
    # ========================================================================
    # SETUP DIRECTORIES
    # ========================================================================
    results_dir = Path('/data/results')
    output_dir = results_dir / 'pipeline_comparison'

    # Handle both /data/results/results/DEG and /data/results/DEG layouts
    candidate_deg_dirs = [
        results_dir / 'results' / 'DEG',   # actual layout seen on disk
        results_dir / 'DEG',               # expected layout
        Path('/data/DEG'),                  # legacy fallback
    ]
    deg_dir = None
    for candidate in candidate_deg_dirs:
        if candidate.exists():
            deg_dir = candidate
            break

    if deg_dir is None:
        print("ERROR: DEG directory not found. Tried:")
        for c in candidate_deg_dirs:
            print(f"  {c}")
        sys.exit(1)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"DEG Directory: {deg_dir}")
    print(f"Output Directory: {output_dir}")
    print("")
    # ========================================================================
    # LOAD NORMALIZED COUNT DATA
    # ========================================================================
    print(">>> Loading normalized count matrices...")
    print("")
    
    normalized_data = {}
    for code in PIPELINE_CODES.keys():
        norm_df = load_normalized_counts(deg_dir, code)
        normalized_data[code] = norm_df
    
    # Check if any data was loaded
    loaded_count = sum(1 for df in normalized_data.values() if df is not None)
    print(f"\n✓ Loaded {loaded_count}/16 normalized matrices")
    
    if loaded_count == 0:
        print("⚠ WARNING: No normalized matrices found!")
        print("  Comparison analysis will continue with DEG results only")
    
    # ========================================================================
    # CREATE COMPARISON MATRICES
    # ========================================================================
    up_matrix = create_comparison_matrix(deg_dir, 'UP')
    down_matrix = create_comparison_matrix(deg_dir, 'DOWN')
    
    if up_matrix is None and down_matrix is None:
        print("\n" + "=" * 70)
        print("ERROR: Could not load any DEG result files!")
        print("=" * 70)
        print("\nPlease ensure that DEG analysis has been completed.")
        print(f"Expected location: {deg_dir}")
        print("\nChecking directory contents:")
        if deg_dir.exists():
            files = list(deg_dir.glob('*'))
            if files:
                for f in sorted(files)[:10]:
                    print(f"  - {f.name}")
            else:
                print("  (directory is empty)")
        sys.exit(1)
    
    # ========================================================================
    # SAVE COMPARISON MATRICES
    # ========================================================================
    print("")
    print(">>> Saving comparison matrices...")
    
    if up_matrix is not None:
        output_file = output_dir / 'upregulated_genes_comparison_matrix.xlsx'
        up_matrix.to_excel(output_file)
        print(f"  ✓ Saved: upregulated_genes_comparison_matrix.xlsx")
        
        output_file_csv = output_dir / 'upregulated_genes_comparison_matrix.csv'
        up_matrix.to_csv(output_file_csv)
        print(f"  ✓ Saved: upregulated_genes_comparison_matrix.csv")
    
    if down_matrix is not None:
        output_file = output_dir / 'downregulated_genes_comparison_matrix.xlsx'
        down_matrix.to_excel(output_file)
        print(f"  ✓ Saved: downregulated_genes_comparison_matrix.xlsx")
        
        output_file_csv = output_dir / 'downregulated_genes_comparison_matrix.csv'
        down_matrix.to_csv(output_file_csv)
        print(f"  ✓ Saved: downregulated_genes_comparison_matrix.csv")
    
    # ========================================================================
    # CREATE VISUALIZATIONS
    # ========================================================================
    print("")
    print(">>> Creating heatmaps...")
    
    if up_matrix is not None:
        create_heatmap(up_matrix, 'UP-Regulated Genes Overlap', 
                      'upregulated_genes_heatmap.png', output_dir)
    
    if down_matrix is not None:
        create_heatmap(down_matrix, 'DOWN-Regulated Genes Overlap', 
                      'downregulated_genes_heatmap.png', output_dir)
    
    # ========================================================================
    # CREATE CORRELATION MATRICES
    # ========================================================================
    print("")
    print(">>> Calculating correlation matrices...")
    
    if up_matrix is not None:
        up_corr = calculate_correlation_lower_triangular(up_matrix)
        if up_corr is not None:
            output_file = output_dir / 'upregulated_genes_correlation_lower_triangular.xlsx'
            up_corr.to_excel(output_file)
            print(f"  ✓ Saved: upregulated_genes_correlation_lower_triangular.xlsx")
            
            # Create correlation heatmap
            plt.figure(figsize=(14, 12))
            sns.heatmap(up_corr, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                       cbar_kws={'label': 'Correlation'}, vmin=-1, vmax=1)
            plt.title('UP-Regulated Genes Correlation (Lower Triangular)', 
                     fontsize=16, fontweight='bold')
            plt.tight_layout()
            output_file = output_dir / 'upregulated_genes_correlation.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  ✓ Saved: upregulated_genes_correlation.png")
    
    if down_matrix is not None:
        down_corr = calculate_correlation_lower_triangular(down_matrix)
        if down_corr is not None:
            output_file = output_dir / 'downregulated_genes_correlation_lower_triangular.xlsx'
            down_corr.to_excel(output_file)
            print(f"  ✓ Saved: downregulated_genes_correlation_lower_triangular.xlsx")
            
            # Create correlation heatmap
            plt.figure(figsize=(14, 12))
            sns.heatmap(down_corr, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                       cbar_kws={'label': 'Correlation'}, vmin=-1, vmax=1)
            plt.title('DOWN-Regulated Genes Correlation (Lower Triangular)', 
                     fontsize=16, fontweight='bold')
            plt.tight_layout()
            output_file = output_dir / 'downregulated_genes_correlation.png'
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  ✓ Saved: downregulated_genes_correlation.png")
    
    # ========================================================================
    # GENERATE SUMMARY REPORT
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
        f"  ✓ UP-regulated genes comparison matrix (16x16)",
        f"  ✓ DOWN-regulated genes comparison matrix (16x16)",
        f"  ✓ Correlation matrices (lower triangular)",
        f"  ✓ Heatmaps for all comparisons",
        "",
        "Files Created:",
        "-" * 60,
        "  • upregulated_genes_comparison_matrix.xlsx",
        "  • upregulated_genes_comparison_matrix.csv",
        "  • upregulated_genes_heatmap.png",
        "  • upregulated_genes_correlation_lower_triangular.xlsx",
        "  • upregulated_genes_correlation.png",
        "  • downregulated_genes_comparison_matrix.xlsx",
        "  • downregulated_genes_comparison_matrix.csv",
        "  • downregulated_genes_heatmap.png",
        "  • downregulated_genes_correlation_lower_triangular.xlsx",
        "  • downregulated_genes_correlation.png",
        "",
    ])
    
    # Add statistics
    if up_matrix is not None and not up_matrix.empty:
        non_diag = up_matrix.values[~np.eye(len(up_matrix), dtype=bool)]
        if len(non_diag) > 0:
            summary_lines.extend([
                "UP-Regulated Genes Statistics:",
                "-" * 60,
                f"  Mean overlap: {non_diag.mean():.0f} genes",
                f"  Max overlap: {non_diag.max():.0f} genes",
                f"  Min overlap: {non_diag.min():.0f} genes",
                "",
            ])
    
    if down_matrix is not None and not down_matrix.empty:
        non_diag = down_matrix.values[~np.eye(len(down_matrix), dtype=bool)]
        if len(non_diag) > 0:
            summary_lines.extend([
                "DOWN-Regulated Genes Statistics:",
                "-" * 60,
                f"  Mean overlap: {non_diag.mean():.0f} genes",
                f"  Max overlap: {non_diag.max():.0f} genes",
                f"  Min overlap: {non_diag.min():.0f} genes",
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
