# ===============================
# RNA-seq Differential Expression Analysis (GSE49712 ENCODE)
# ===============================

# ---- . countsقراءة بيانات الـ  ----

EXprSet <- read.delim(
  "C:/Users/HP/Downloads/GSE49712_ENCODE_HTSeq.txt/GSE49712_ENCODE_HTSeq.txt")
# تحويلها لمصفوفة للتعامل مع DESeq2

CountData <- as.matrix(EXprSet)
# ---- 2. إعداد معلومات العينات ----
ColData <-data.frame(
  condition = factor(c(
    rep("GM12892" ,3),
    rep("H1.hES" ,4),
    rep("MCF.7" ,3)
  )),
  row.names = colnames(CountData))
# ----  إنشاء DESeq2 dataset ----
library(DESeq2)
dds = DESeqDataSetFromMatrix(
  countData = CountData,
  colData = ColData,
  design = ~condition
)
dds =DESeq(dds)
res = results(dds)
res =res[order(res$padj),]
sig.genes = res[!is.na(res$padj)& res$padj < .05 , ]
write.csv(
  as.data.frame(sig.genes),
  "C:/Users/HP/Documents/RNAseq_Project/finally.csv"
)

# رسم الـ Volcano plot


library(ggplot2)
volcano = as.data.frame(res)
volcano$significant =ifelse(!is.na(volcano$padj) & is.finite(volcano$padj) & (volcano$padj)< .05 ,"yes" ,"no")
volcano$pointsize = -log10(volcano$padj)
volcano$pointsize = pmax(pmin(volcano$pointsize ,3),1)

ggplot(volcano, aes( x = log2FoldChange, y=-log10(padj) ,colour = significant))+
  geom_point(aes(size= pointsize) ,alpha = .6)+
  scale_color_manual(values=c("no"= "grey" ,"yes"= "red"))+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust= 0.5, face = "bold", size = 10 ,colour = "black"),
    axis.title = element_text(face = "bold" , size = 12)
  )+
  labs(title = "Differential Gene Expression Analysis of GSE49712 (ENCODE RNA-seq)",
       x="log2 Fold Change" ,
       y= "-log10 adjusted p-value")

sig.genes = sig.genes[order(abs(sig.genes$log2FoldChange),decreasing = T),]
top.gene =head(sig.genes,500)

# ---- . استخراج القراءات المعنوية للـ Heatmap ----
normcounts = counts(dds,normalized =T)
sigcounts = normcounts[row.names(normcounts) %in% row.names(top.gene), , drop = FALSE]
sigcounts =sigcounts[rowSums(sigcounts)>0,]

# ----  تقيس الجينات لكل صف(z-score) ----
sigcounts_scaled = t(scale(t(sigcounts)))
sigcounts_scaled= sigcounts_scaled[!apply(is.na(sigcounts_scaled)|is.infinite(sigcounts_scaled),1,any),]


# رسم الـ  heatmap
library(pheatmap)
pheatmap(sigcounts_scaled,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = T,color = colorRampPalette(c("blue" ,"white","red"))(100),
         fontsize = 10,
         main = "Differential Gene Expression Analysis of GSE49712(ENCODE RNA-seq)")


