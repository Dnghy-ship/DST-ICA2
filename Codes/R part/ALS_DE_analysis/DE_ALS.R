#参考数据集
secreted <- read.table(
  "Secreted_protein_ref_ENSEMBLE_intersection.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)[,1]

lr <- read.csv("Secreted_to_Receptors_Ramilowski.csv")
ligands   <- unique(lr$ligand_ENSG)
receptors <- unique(lr$receptor_ENSG)
lr_genes <- unique(c(ligands, receptors))

#读入实验数据
rna <- read.csv("ALS_RNA_GSE287256_spinal_cord.csv", row.names = 1)

#clean ENSG
clean_ensg <- function(x){
  x <- as.character(x)
  x <- gsub("\\..*", "", x)   # 去掉 ENSG 的版本号
  x <- gsub("\\|.*", "", x)   # 去掉 symbol 等
  return(x)
}

rownames(rna) <- clean_ensg(rownames(rna))

#筛选
expr_secreted <- rna[rownames(rna) %in% secreted, ]
expr_receptor <- rna[rownames(rna) %in% receptors, ]
expr_lr <- rna[rownames(rna) %in% lr_genes, ]

#输出保存
write.csv(expr_secreted, "ALS_secreted_expression.csv")
write.csv(expr_receptor, "ALS_receptor_expression.csv")
write.csv(expr_lr, "ALS_LR_expression.csv")

#准备差异分析
#检查
head(expr_secreted[,1:4])

#读入
expr <- read.csv("ALS_secreted_expression.csv", row.names = 1)

#构建metadata表
samples <- colnames(expr)

condition <- ifelse(grepl("^CT", samples), "Control", "ALS")

metadata <- data.frame(
  row.names = samples,
  condition = factor(condition, levels = c("Control", "ALS"))
)
metadata

#DESeq2
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = round(expr),
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "ALS", "Control"))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df[complete.cases(res_df), ]

#筛选差异基因
sig <- res_df[
  abs(res_df$log2FoldChange) >= 1 & res_df$padj < 0.05,
]
up_secreted   <- sig[sig$log2FoldChange > 1, ]
down_secreted <- sig[sig$log2FoldChange < -1, ]

#火山图
library(ggplot2)

res_df$diff <- "Not Sig"
res_df$diff[res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "Up"
res_df$diff[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = diff)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(
    title = "ALS vs Control — Secreted Genes",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

#Top 30 genes
library(pheatmap)

top_genes <- head(sig[order(sig$padj), "gene"], 30)

heatmap_mat <- expr[top_genes, ]

pheatmap(
  log2(heatmap_mat + 1),
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  annotation_col = metadata,
  main = "Top 30 Differential Secreted Genes (ALS)"
)

#阶段结果导出
write.csv(res_df, "ALS_secreted_DE_full.csv", row.names = FALSE)
write.csv(sig, "ALS_secreted_DE_sig.csv", row.names = FALSE)
write.csv(up_secreted, "ALS_secreted_up.csv", row.names = FALSE)
write.csv(down_secreted, "ALS_secreted_down.csv", row.names = FALSE)

