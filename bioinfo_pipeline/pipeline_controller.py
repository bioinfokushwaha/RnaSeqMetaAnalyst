#!/usr/bin/env python3

import argparse
import subprocess
import os
import csv
from datetime import datetime

def run_command(label, command, log_dir, summary):
    log_file = os.path.join(log_dir, f"{label}.log")
    with open(log_file, 'w') as log:
        try:
            subprocess.run(command, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
            summary.append([label, "SUCCESS", log_file])
        except subprocess.CalledProcessError:
            summary.append([label, "FAILED", log_file])
            print(f"[ERROR] Step '{label}' failed. Check {log_file}")

def main():
    parser = argparse.ArgumentParser(description="Unified Bioinformatics Pipeline Controller")
    parser.add_argument("--aligners", type=str, help="Comma-separated: hisat2,star,bowtie2")
    parser.add_argument("--quantifiers", type=str, help="Comma-separated: htseq,featurecounts,rsem")
    parser.add_argument("--deanalysis", type=str, help="Comma-separated: deseq2,edger")
    parser.add_argument("--mode", choices=["PE", "SE"], default="PE", help="Sequencing mode (Paired/Single End)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--ref", required=True, help="Path to genome reference folder")
    parser.add_argument("--reads", required=True, help="Path to read files folder")
    parser.add_argument("--outdir", default="output", help="Main output directory")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    log_dir = os.path.join(args.outdir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    summary = []
    aligners = args.aligners.split(",") if args.aligners else []
    quantifiers = args.quantifiers.split(",") if args.quantifiers else []
    de_analysis = args.deanalysis.split(",") if args.deanalysis else []

    print("=== Starting Bioinformatics Pipeline ===")
    print(f"Aligners: {aligners}")
    print(f"Quantifiers: {quantifiers}")
    print(f"Differential Expression: {de_analysis}")
    print(f"Mode: {args.mode}, Threads: {args.threads}")
    print(f"Reference: {args.ref}, Reads: {args.reads}")
    print(f"Output Directory: {args.outdir}")

    #################### ALIGNMENT ####################
    for aligner in aligners:
        script = f"scripts/run_pipeline_{aligner}.sh"
        if os.path.exists(script):
            cmd = f"bash {script} {args.ref} {args.reads} {args.mode} {args.threads}"
            run_command(f"{aligner}_align", cmd, log_dir, summary)
        else:
            summary.append([f"{aligner}_align", "SCRIPT_NOT_FOUND", script])

    #################### QUANTIFICATION ####################
    for quant in quantifiers:
        script = f"scripts/run_{quant}.sh"
        if os.path.exists(script):
            cmd = f"bash {script} {args.ref} {args.mode} {args.threads}"
            run_command(f"{quant}_quant", cmd, log_dir, summary)
        else:
            summary.append([f"{quant}_quant", "SCRIPT_NOT_FOUND", script])

    #################### DIFFERENTIAL EXPRESSION ####################
    os.makedirs(os.path.join(args.outdir, "DE_Analysis"), exist_ok=True)

    for de in de_analysis:
        if de == "deseq2":
            cmd = ("Rscript scripts/de_deseq2.R Sampleinfo.txt "
                   "Quantification/featurecounts/FC_Count.txt "
                   f"{args.outdir}/DE_Analysis/DESeq2_logFC.txt {args.outdir}/DE_Analysis/DESeq2_Pval.txt")
            run_command("deseq2", cmd, log_dir, summary)

        elif de == "edger":
            cmd = ("Rscript scripts/de_edger.R Sampleinfo.txt "
                   "Quantification/featurecounts/FC_Count.txt "
                   f"{args.outdir}/DE_Analysis/edgeR_logFC.txt {args.outdir}/DE_Analysis/edgeR_Pval.txt")
            run_command("edger", cmd, log_dir, summary)

        else:
            summary.append([de, "UNKNOWN_DE_METHOD", ""])

    #################### SUMMARY OUTPUT ####################
    summary_file = os.path.join(args.outdir, "pipeline_summary.csv")
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Step", "Status", "LogFile"])
        writer.writerows(summary)

    print(f"\n=== Pipeline completed ===\nSummary written to: {summary_file}")

if __name__ == "__main__":
    main()

