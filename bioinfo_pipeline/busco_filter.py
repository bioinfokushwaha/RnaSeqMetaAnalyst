#!/usr/bin/env python3
# ===================================================================
# busco_filter.py — Filter GFF3 / GTF to BUSCO core genes
#                   + generate housekeeping_genes.xlsx
# ===================================================================
# Environment : rnaseq-busco-py  (Python 3.11, gffpandas, pandas 2.x)
#
# Steps:
#   1.  Load BUSCO complete single-copy protein IDs
#   2.  Read GFF3
#   3.  Extract protein_id from GFF3 attributes
#   4.  Filter GFF3 rows to BUSCO matches
#   5.  Read GTF
#   6.  Parse gene_id and protein_id from GTF
#       NOTE: handles NCBI doubled-quote format  gene_id ""KCNE2""
#             as well as standard format         gene_id "KCNE2"
#   7.  Filter GTF rows whose protein_id is in BUSCO list
#   8.  Merge GFF3 + GTF into summary table
#   9.  Save busco_filtered.gff, busco_filtered.gtf, busco_matched.tsv
#   10. Generate housekeeping_genes.xlsx
#       - extracts gene_id symbols (e.g. KCNE2, DONSON) from GTF
#       - matches via protein_id to confirm they are BUSCO genes
#       - header column: gene_id  (matches R scripts hkg_data[[1]])
#       - saves to <out-dir>/housekeeping_genes.xlsx
#       - auto-copies to /data/housekeeping_genes.xlsx so DESeq2
#         and edgeR pick it up without any manual step
#
# Usage inside container:
#   micromamba run -n rnaseq-busco-py \
#     python /opt/project/scripts/busco_filter.py \
#       --gff      /data/genome/genomic.gff \
#       --gtf      /data/genome/genomic.gtf \
#       --busco    /data/busco_output/busco_core_genes1.txt \
#       --out-dir  /data/busco_output/filtered \
#       --data-dir /data
# ===================================================================

import argparse
import os
import re
import shutil
import sys

import gffpandas.gffpandas as gffpd
import pandas as pd


# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter GFF3 and GTF to BUSCO complete single-copy genes "
            "and generate housekeeping_genes.xlsx for DESeq2/edgeR."
        )
    )
    parser.add_argument("--gff",      default="genomic.gff",
                        help="Input GFF3 file")
    parser.add_argument("--gtf",      default="genomic.gtf",
                        help="Input GTF file")
    parser.add_argument("--busco",    default="busco_core_genes1.txt",
                        help="BUSCO complete gene list — TSV, protein IDs in column 0")
    parser.add_argument("--out-dir",  dest="out_dir",  default=".",
                        help="Output directory for all result files")
    parser.add_argument("--data-dir", dest="data_dir", default="/data",
                        help="Container /data root where housekeeping_genes.xlsx "
                             "is auto-copied for DESeq2/edgeR (default: /data)")
    return parser.parse_args()


# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
def extract_protein_gff(attr: str) -> str:
    """Extract protein_id values from a GFF3 attributes string."""
    attr = str(attr)
    proteins = []
    # Standard GFF3: protein_id=NP_001070552.1
    proteins.extend(re.findall(r'protein_id=([^;,\s]+)', attr))
    # Genbank tag: Genbank:XP_... or NP_...
    proteins.extend(re.findall(r'Genbank:(XP_[^,;\s]+|NP_[^,;\s]+)', attr))
    return ",".join(p.strip() for p in proteins if p.strip())


def parse_gtf_attr(attr: str, key: str) -> str:
    """
    Extract a key value from a GTF attributes string.
    Handles two formats:
      NCBI doubled-quote:  gene_id ""KCNE2"";
      Standard:            gene_id "KCNE2";
    """
    attr = str(attr)
    # NCBI doubled-quote format (e.g. from BestRefSeq GTFs)
    m = re.search(rf'{key}\s+""([^"]+)""', attr)
    if m:
        return m.group(1).strip()
    # Standard GTF format
    m = re.search(rf'{key}\s+"([^"]+)"', attr)
    return m.group(1).strip() if m else ""


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    output_table  = os.path.join(args.out_dir, "busco_matched.tsv")
    output_gff    = os.path.join(args.out_dir, "busco_filtered.gff")
    output_gtf    = os.path.join(args.out_dir, "busco_filtered.gtf")
    output_hkg    = os.path.join(args.out_dir, "housekeeping_genes.xlsx")
    data_hkg_copy = os.path.join(args.data_dir, "housekeeping_genes.xlsx")

    # ------------------------------------------------------------------
    # STEP 1: Load BUSCO protein IDs
    # ------------------------------------------------------------------
    print(">>> Step 1: Loading BUSCO complete gene protein IDs ...")
    if not os.path.exists(args.busco):
        sys.exit(f"ERROR: BUSCO file not found: {args.busco}")

    busco_df  = pd.read_csv(args.busco, sep="\t", header=None)
    busco_ids = set(busco_df[0].astype(str).str.strip())
    print(f"  Loaded {len(busco_ids):,} protein IDs")

    # ------------------------------------------------------------------
    # STEP 2: Read GFF3
    # ------------------------------------------------------------------
    print(">>> Step 2: Reading GFF3 ...")
    if not os.path.exists(args.gff):
        sys.exit(f"ERROR: GFF file not found: {args.gff}")

    annotation = gffpd.read_gff3(args.gff)
    gff_df     = annotation.df.copy()
    print(f"  GFF3 rows loaded: {len(gff_df):,}")

    # ------------------------------------------------------------------
    # STEP 3: Extract protein IDs from GFF3 attributes
    # ------------------------------------------------------------------
    print(">>> Step 3: Extracting protein IDs from GFF3 ...")
    gff_df["protein_id"]   = gff_df["attributes"].apply(extract_protein_gff)
    gff_df["protein_list"] = gff_df["protein_id"].str.split(",")
    gff_df = gff_df.explode("protein_list")
    gff_df["protein_list"] = gff_df["protein_list"].str.strip()

    # ------------------------------------------------------------------
    # STEP 4: Filter GFF3 to BUSCO matches
    # ------------------------------------------------------------------
    print(">>> Step 4: Filtering GFF3 to BUSCO protein IDs ...")
    gff_busco = gff_df[gff_df["protein_list"].isin(busco_ids)].copy()
    print(f"  GFF3 rows matched: {len(gff_busco):,}")
    print(f"  Unique proteins  : {gff_busco['protein_list'].nunique():,}")

    # ------------------------------------------------------------------
    # STEP 5: Read GTF
    # ------------------------------------------------------------------
    print(">>> Step 5: Reading GTF ...")
    if not os.path.exists(args.gtf):
        sys.exit(f"ERROR: GTF file not found: {args.gtf}")

    gtf = pd.read_csv(
        args.gtf,
        sep="\t",
        comment="#",
        header=None,
        names=["seqid","source","type","start","end",
               "score","strand","phase","attributes"]
    )
    print(f"  GTF rows loaded: {len(gtf):,}")

    # ------------------------------------------------------------------
    # STEP 6: Parse gene_id and protein_id from GTF
    #         Uses universal regex for doubled-quote NCBI format
    # ------------------------------------------------------------------
    print(">>> Step 6: Parsing gene_id and protein_id from GTF ...")
    print("  (handles NCBI doubled-quote format: gene_id \"\"KCNE2\"\")")

    gtf["protein_id"] = gtf["attributes"].apply(
        lambda x: parse_gtf_attr(x, "protein_id")
    )
    gtf["gene_id_sym"] = gtf["attributes"].apply(
        lambda x: parse_gtf_attr(x, "gene_id")
    )

    # Quick sanity check
    sample_mask = gtf["protein_id"].str.strip() != ""
    n_with_protein = sample_mask.sum()
    print(f"  GTF rows with protein_id : {n_with_protein:,}")

    sample_mask2 = gtf["gene_id_sym"].str.strip() != ""
    n_with_gene = sample_mask2.sum()
    print(f"  GTF rows with gene_id    : {n_with_gene:,}")

    if n_with_protein > 0:
        example = gtf.loc[sample_mask, ["protein_id","gene_id_sym"]].head(3)
        print(f"  Sample parse:")
        for _, row in example.iterrows():
            print(f"    protein_id={row['protein_id']}  gene_id={row['gene_id_sym']}")

    # ------------------------------------------------------------------
    # STEP 7: Filter GTF rows whose protein_id is in BUSCO list
    # ------------------------------------------------------------------
    print(">>> Step 7: Filtering GTF to BUSCO protein IDs ...")
    gtf_busco = gtf[gtf["protein_id"].isin(busco_ids)].copy()
    print(f"  GTF rows matched : {len(gtf_busco):,}")
    print(f"  Unique proteins  : {gtf_busco['protein_id'].nunique():,}")
    print(f"  Unique gene IDs  : {gtf_busco['gene_id_sym'].nunique():,}")

    # ------------------------------------------------------------------
    # STEP 8: Merge GFF3 + GTF → summary table
    # ------------------------------------------------------------------
    print(">>> Step 8: Merging GFF3 and GTF annotations ...")
    merged = pd.merge(
        gff_busco,
        gtf_busco[["protein_id", "gene_id_sym"]].drop_duplicates(),
        left_on="protein_list",
        right_on="protein_id",
        how="inner"
    )
    print(f"  Merged rows: {len(merged):,}")

    # ------------------------------------------------------------------
    # STEP 9: Save filtered GFF3, GTF, TSV
    # ------------------------------------------------------------------
    print(">>> Step 9: Saving filtered annotations ...")

    # TSV summary table
    merged.to_csv(output_table, sep="\t", index=False)
    print(f"  ✓ TSV table      : {output_table}")

    # Filtered GFF3
    gff_out    = gff_df[gff_df["protein_list"].isin(busco_ids)].copy()
    clean_cols = [c for c in gff_out.columns
                  if c not in ("protein_id", "protein_list")]
    annotation.df = gff_out[clean_cols].reset_index(drop=True)
    annotation.to_gff3(output_gff)
    print(f"  ✓ Filtered GFF3  : {output_gff}")

    # Filtered GTF (drop helper columns)
    gtf_busco.drop(
        columns=["protein_id", "gene_id_sym"], errors="ignore"
    ).to_csv(output_gtf, sep="\t", index=False)
    print(f"  ✓ Filtered GTF   : {output_gtf}")

    # ------------------------------------------------------------------
    # STEP 10: Generate housekeeping_genes.xlsx
    #
    # Strategy:
    #   - Take each unique protein_id that matched BUSCO
    #   - Look up its gene_id symbol from the GTF (e.g. DONSON, KCNE2)
    #   - Result: one row per unique protein_id → gene_id pair
    #   - Column header: gene_id  (matches R scripts hkg_data[[1]])
    #   - If no gene_id symbol found → keep protein_id as fallback
    # ------------------------------------------------------------------
    print(">>> Step 10: Generating housekeeping_genes.xlsx ...")

    # Build protein_id → gene_id mapping from filtered GTF
    # One protein can appear many times (one row per exon/CDS)
    # Take the first non-empty gene_id for each protein_id
    protein_to_gene = (
        gtf_busco[["protein_id", "gene_id_sym"]]
        .pipe(lambda df: df[df["protein_id"].str.strip() != ""])
        .drop_duplicates(subset="protein_id", keep="first")
        .set_index("protein_id")["gene_id_sym"]
        .to_dict()
    )

    # Build the HKG table: iterate over BUSCO protein IDs
    rows = []
    no_symbol = []
    for pid in sorted(busco_ids):
        symbol = protein_to_gene.get(pid, "").strip()
        if symbol:
            rows.append({"gene_id": symbol, "protein_id": pid})
        else:
            # No GTF match — use protein ID itself as fallback
            rows.append({"gene_id": pid, "protein_id": pid})
            no_symbol.append(pid)

    hkg_df = (
        pd.DataFrame(rows)
        .drop_duplicates(subset="gene_id")
        .sort_values("gene_id")
        .reset_index(drop=True)
    )

    print(f"  BUSCO proteins total     : {len(busco_ids):,}")
    print(f"  Mapped to gene symbols   : {len(busco_ids) - len(no_symbol):,}")
    print(f"  Fallback (protein ID)    : {len(no_symbol):,}")
    print(f"  Total HKG entries        : {len(hkg_df):,}")

    if len(hkg_df) > 0:
        print(f"  Sample entries:")
        for _, row in hkg_df.head(5).iterrows():
            print(f"    gene_id={row['gene_id']}  protein_id={row['protein_id']}")

    # Save Excel — only gene_id column (what R scripts read)
    # protein_id kept as second column for traceability
    hkg_df.to_excel(output_hkg, index=False, header=False)
    print(f"  ✓ Saved to out-dir       : {output_hkg}")

    # Auto-copy to /data/ so DESeq2 and edgeR find it automatically
    if os.path.isdir(args.data_dir):
        shutil.copy2(output_hkg, data_hkg_copy)
        print(f"  ✓ Auto-copied to /data/  : {data_hkg_copy}")
    else:
        print(f"  ⚠ data-dir not found ({args.data_dir})")
        print(f"    Manually copy {output_hkg} to your results directory")

    # ------------------------------------------------------------------
    # FINAL SUMMARY
    # ------------------------------------------------------------------
    print("")
    print("=" * 55)
    print("BUSCO FILTER COMPLETE")
    print("=" * 55)
    print(f"  BUSCO proteins loaded      : {len(busco_ids):,}")
    print(f"  GFF3 rows matched          : {len(gff_busco):,}")
    print(f"  GTF rows matched           : {len(gtf_busco):,}")
    print(f"  Unique proteins in GFF3    : {gff_busco['protein_list'].nunique():,}")
    print(f"  Unique proteins in GTF     : {gtf_busco['protein_id'].nunique():,}")
    print(f"  Gene symbols resolved      : {len(busco_ids) - len(no_symbol):,}")
    print(f"  Housekeeping genes written : {len(hkg_df):,}")
    print("")
    print("Output files:")
    print(f"  {output_table}")
    print(f"  {output_gff}")
    print(f"  {output_gtf}")
    print(f"  {output_hkg}")
    print(f"  {data_hkg_copy}  ← auto-picked up by DESeq2 & edgeR")
    print("=" * 55)


if __name__ == "__main__":
    main()
