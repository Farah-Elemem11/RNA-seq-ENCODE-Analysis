# RNA-Seq Differential Expression Analysis | MCF.7 vs GM12892

## Project Overview

This project demonstrates a complete RNA-seq differential expression analysis workflow using the DESeq2 package in R. The analysis compares the breast cancer cell line **MCF.7** with the normal lymphoblastoid cell line **GM12892** to identify significantly differentially expressed genes (DEGs). The **H1.hESC** embryonic stem cell line was included for exploratory visualization in the heatmap.

The workflow includes data preprocessing, normalization, statistical testing, differential expression analysis, and publication-quality visualizations to reveal transcriptomic differences between the analyzed cell lines.

## Tools & Technologies

**Programming Language**

* R

**Core Packages**

* DESeq2 – Differential expression analysis
* ggplot2 – Volcano plot visualization
* pheatmap – Clustered heatmap visualization
* RColorBrewer – Color palette optimization

## Workflow Pipeline

1. Loaded raw HTSeq count data.
2. Constructed sample metadata and defined the experimental design.
3. Filtered low-abundance genes (counts ≥ 10).
4. Normalized gene counts and performed differential expression analysis using DESeq2.
5. Identified significant DEGs using:

   * Adjusted *P*-value < 0.05
   * |log2 Fold Change| > 1
6. Generated a volcano plot to visualize significantly regulated genes.
7. Created a Z-score heatmap of the top 50 differentially expressed genes.

## Key Results

* Identified significantly upregulated and downregulated genes between **MCF.7** and **GM12892**.
* Generated publication-quality volcano plots highlighting statistically significant DEGs.
* Produced normalized heatmaps demonstrating expression patterns of the top 50 DEGs.
* Observed clear clustering of samples, supporting the biological consistency of the dataset.

## Repository Structure

* **GSE49712_ENCODE_HTSeq.txt** – Raw HTSeq count matrix.
* **RNAseq_Final_Analysis.R** – Complete analysis workflow.
* **RNAseq_Final_Results.csv** – Differential expression results.
* **Plots/** – Volcano plot and heatmap figures.
