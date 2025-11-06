RNA-seq Differential Expression Analysis (GSE49712 - ENCODE)

This project performs RNA-seq differential gene expression analysis using the DESeq2 package in R, applied to the publicly available dataset GSE49712 (ENCODE project).

Overview
The aim of this analysis is to identify genes that are differentially expressed among three human cell lines:

* GM12892
* H1.hES
* MCF.7

Workflow Steps

1. Import raw count data from HTSeq files
2. Construct the DESeq2 dataset
3. Perform normalization and differential expression analysis
4. Extract significantly expressed genes (adjusted p < 0.05)
5. Visualize results using Volcano Plot and Heatmap


Tools & Libraries Used

* R / RStudio
* DESeq2
* ggplot2
* pheatmap



Project Structure

RNAseq_DESeq2_GSE49712/
│
├── README.md                (Project description - this file)
│
├── data/                    (Raw data files)
│   └── GSE49712_ENCODE_HTSeq.txt
│
├── scripts/                 (Analysis scripts)
│   └── DESeq2_analysis.R
│
├── results/                 (CSV results and outputs)
│   └── finally.csv
│
└── figures/                 (Generated figures)
├── volcano_plot.png
└── heatmap.png



Key Outputs

* finally.csv → List of significantly differentially expressed genes
* volcano_plot.png → Visual summary of gene significance vs fold change
* heatmap.png → Expression pattern of top 500 significant genes

Author
Farah Mohamed Elimam Mohamed
M.Sc. Preliminary Student in Biochemistry, Mansoura University
Egypt


Notes
This project represents a small part of my practical training in RNA-seq data analysis and bioinformatics using R.
Further steps will include functional enrichment and pathway analysis in upcoming projects.





