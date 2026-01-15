Author: Farah Mohamed Elemam > Role:Biochemist & Bioinformatics Researcher 
 
 RNA-Seq Differential Expression Analysis | MCF.7 vs GM12892

 Project Overview

This project performs a comprehensive Transcriptomic Analysis to identify Differentially Expressed Genes (DEGs) between different cell lines: MCF.7 (Breast Cancer), GM12892 (Lymphoblastoid - Normal Control), and H1.hESC (Embryonic Stem Cells).

The analysis aims to uncover the molecular signatures that distinguish cancerous cells from normal and stem cells using the DESeq2 framework.

 Tools & Technologies

Language: R
Core Packages: DESeq2: For differential expression testing.
ggplot2: For high-quality data visualization (Volcano Plots).
pheatmap: For clustered gene expression heatmaps.
RColorBrewer: For optimized color palettes.


 Workflow Pipeline

1. Data Pre-processing: Loading raw HTSeq counts and filtering out low-abundance genes (counts < 10).
2.Metadata Construction: Defining experimental conditions and setting the reference baseline (GM12892).
3.Statistical Analysis: Normalization and fitting the Negative Binomial distribution using DESeq2.
4. Significance Filtering: Identifying DEGs based on:
 Adjusted P-value < 0.05
 |log2FoldChange| > 1


6. Visualization: Generating Volcano plots to show the global distribution of DEGs and Heatmaps to show the Z-score scaling of the top 50 genes.

  Key Results

 Up-regulated Genes: Identification of genes overexpressed in the MCF.7 cancer line, potentially acting as oncogenic drivers.
 Down-regulated Genes: Identification of genes suppressed in cancer, often associated with normal cellular functions or tumor suppression.
 Clustering: Clear separation between cell types in the heatmap, validating the biological consistency of the samples.

  Repository Structure

GSE49712_ENCODE_HTSeq.txt: The raw input count matrix.
RNAseq_Final_Analysis.R: The complete, optimized R script.
RNAseq_Final_Results.csv: Output table containing log2FoldChange and P-values.
Plots: Directory containing the Volcano Plot and Heatmap images.


onnect with Me
LinkedIn: https://www.linkedin.com/in/farah-elemam-107969323



