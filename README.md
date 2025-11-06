# 🧬 RNA-seq Differential Expression Analysis (GSE49712 - ENCODE)

This project performs **RNA-seq differential gene expression analysis** using the **DESeq2** package in R, applied on the publicly available dataset **GSE49712 (ENCODE project)**.

---

 📘 Overview
The aim of this analysis is to identify genes that are **differentially expressed** among three human cell lines:
- GM12892
- H1.hES
- MCF.7**

The workflow includes:
1. Importing raw count data from HTSeq files  
2. Constructing the DESeq2 dataset  
3. Performing normalization and differential expression analysis  
4. Extracting significantly expressed genes (adjusted *p* < 0.05)  
5. Visualization using Volcano Plot and Heatmap


 Tools & Libraries Used
- R / RStudio
- DESeq2
- ggplot2
- pheatmap



## 📂 Project Structure
