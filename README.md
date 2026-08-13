# gse326385-dge-and-crispr-design
A complete bioinformatics workflow evaluating a Free Fatty Acid (FFA)-induced metabolic stress model (GSE326385). Includes a raw RNA-Seq count differential gene expression (DGE) pipeline using DESeq2 in R, publication-quality volcano plots, and a targeted CRISPR/Cas9 sgRNA knockout architecture for the upregulated master driver EGR2.

Project Overview
This repository contains a unified bioinformatic workflow executing two core analytical milestones:
1. Differential Gene Expression (DGE) parsing using DESeq2 in R to evaluate a Free Fatty Acid (FFA)-induced metabolic cell stress model.
2. CRISPR/Cas9 Guide RNA Engineering via CHOPCHOP to isolate viable gene knockouts targeting the master regulatory driver EGR2.

Task 1: Differential Gene Expression (DGE)
Methodology
Pipeline Framework: DESeq2 empirical Bayes analysis.
Data Dimension: 18 high-throughput transcriptomic replicates (9 Control baseline vs 9 FFA-Treated).
Significance Thresholds: Adjusted P-value (FDR) < 0.05, |Log2 Fold Change| > 1.0.

Visual Output Distribution
The volcano plot below demonstrates global transcriptomic variations under treatment stress, distinguishing significantly upregulated pathways (Red) and downregulated baselines (Blue). Detailed breakdown of the plot is below:
1.	Upregulated Genes (Red Dots):
•	These genes have a Log₂FC > 1 (meaning their expression more than doubled) and an adjusted p-value < 0.05.

•	Biological Meaning: These genes are actively triggered by FFA treatment. They represent paths involved in cellular stress, fat storage dynamics, inflammation, or lipotoxicity. Your chosen knockout target, EGR2, sits within this group.

2.	Downregulated Genes (Blue Dots):
•	These genes have a Log₂FC < -1 (expression dropped by more than half) and an adjusted p-value < 0.05.

•	Biological Meaning: These are genes that the cell suppresses or shuts down when exposed to fatty acid toxicity.

3.	Not Significant Genes (Grey Dots):
These are genes that either did not change much between the conditions, or had too much variation to be statistically reliable.


Key observations:
•	Highly Dynamic Outlier: That single red dot sitting far out on the top-right corner (around Log₂FC ≈ 5 and -Log₁₀p-adj ≈ 7) shows that gene's expression got increased by over 32-fold (2⁵) under FFA treatment with massive statistical significance. 

•	Asymmetric Distribution: The plot features slightly more upregulated (red) dots than downregulated (blue) dots, showing that FFA treatment acts primarily as a transcriptional activator that forces the liver model cells to change their behaviour.

Task 2: CRISPR/Cas9 Guide RNA Engineering
Design Objective
To systematically silence the master zinc-finger transcription factor EGR2 which was found significantly overexpressed under lipotoxic stress in Task 1.

Guide Design Parameters
Target Software Platform: CHOPCHOP Suite (v3)
Organism Assembly Target: Homo sapiens (hg38)
Enzyme Core Model: Cas9 (PAM sequence requirement: NGG)
Targeting Region Bounds: Coding Exon Intersection Alignment

Selected Guide Candidates

RANK	 Target Sequence (5’ -> 3’)	  PAM	 Genomic Location	 GC Content	  MM0/MM1/MM2/MMI
 1	       TCAAGGTGTCCGGGTCCGAG	    AGG	  chr10:62813207	     65	          0/0/0/0
 2	       GCAAGACGCCGGTGCACGAG	    AGG	  chr10:628132624	     70	          0/0/0/0
