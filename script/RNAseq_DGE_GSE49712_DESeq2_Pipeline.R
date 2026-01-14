# =============================================================
# RNA-seq Differential Expression Analysis (GSE49712 ENCODE)
# Project: Comparison of MCF.7, H1.hESC, and GM12892 Cell Lines
# Author: Farah (Biochemist)
# =============================================================


#--- 1. Loading Essential Libraries ---
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#--- 2. Data Acquisition & Matrix Preparation ---
raw_counts <- read.delim("C:/Users/HP/Downloads/RNAseq_DESeq2_GSE49712/GSE49712_ENCODE_HTSeq.txt", 
                         header = T , check.names = F, sep = "\t", row.names = 1)
count_matrix = as.matrix( raw_counts)

#--- 3. Experimental Design (Metadata) ---
col_Data =data.frame(condition = factor(c( rep("GM12892" ,3) ,rep( "H1.hESC", 4) ,rep("MCV.7" ,3))) ,
row.names = colnames(count_matrix))

#--- 4. DESeq2 Object Construction & Pre-filtering ---
dds = DESeqDataSetFromMatrix(countData =count_matrix , colData  = col_Data , design =  ~ condition  )
keep = rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "GM12892")

#--- 5. Differential Expression Analysis (The Core Pipeline) ---
dds =DESeq(dds)
# Extracting results for the specific comparison (Test vs Reference)
res = results(dds , contrast = c("condition","MCV.7","GM12892"))
res_df = as.data.frame(res)
keep <- rowSums(counts(dds)) >= 10
dds = dds[keep,]

#--- 6. Significance Labeling ---
res_df$status="not_sig"
res_df$status[res_df$padj <0.05 &res_df$log2FoldChange>1] = "UP"
res_df$status[res_df$padj <0.05 &res_df$log2FoldChange< -1] ="DOWN"
res_df=  res_df[order(res_df$padj),]
write.csv(res_df,"RNAseq_Finall_Results.csv")

#--- 7. Volcano Plot Visualization ---
ggplot(res_df, aes(x = log2FoldChange , y =-log10(padj) ,col = status))+
  geom_point(alpha= .4 ,size = 1.5 )+
  scale_colour_manual(values = c("UP" = "red" ,"DOWN" ="blue" ,"not_sig" = "grey"))+
  theme_minimal()+
  geom_vline(xintercept = c(-1,1) ,color ="black" ,linetype = "dashed")+
  geom_hline(yintercept = -log10(.05),color ="black" ,linetype ="dashed")+
  labs(title = "Volcano Plot: MCF.7 vs GM12892",
       x = "log2 Fold Change", y = "-log10 Adjusted P-value")

#--- 8. Heatmap Preparation & Z-score Scaling ---
top_gene_idx=rownames(res_df [1:50,])
normcounts = counts(dds, normalized = T)
sig_counts =normcounts[top_gene_idx,]

#--- 9. Annotated Heatmap Construction ---
sig_count_scaled =t(scale(t(sig_counts)))
sample_info = as.data.frame(colData(dds)[,"condition" , drop =F])
colnames(sample_info) ="Group"

ann_colors= list(Group = c("MCV.7" = "red", "H1.hESC" = "grey", "GM12892" = "blue"))

pheatmap(sig_count_scaled,
         annotation_col = sample_info,
         annotation_colors = ann_colors,
         show_rownames = T,
         show_colnames = F,         
         color = colorRampPalette(c("blue","white", "red"))(50),
         main = "Heatmap of Top 50 DEGs (Z-score)")
  
  
  
  

