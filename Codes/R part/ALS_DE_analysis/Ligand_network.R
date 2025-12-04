#构建L-R网络
#读入
library(dplyr)

deg <- read.csv("ALS_secreted_DE_sig.csv")

head(deg)

lr_ref <- read.csv("Secreted_to_Receptors_Ramilowski.csv")

head(lr_ref)

#筛选
lr_deg <- lr_ref %>%
  inner_join(deg, by = c("ligand_ENSG" = "gene"))

#构建
ALS_LR_network <- lr_deg %>%
  select(
    ligand_ENSG,
    receptor_ENSG,
    log2FoldChange,
    padj
  )
write.csv(ALS_LR_network, "ALS_Ligand_Receptor_Network.csv", row.names = FALSE)

head(ALS_LR_network,20)

#可视化
# 依赖包
library(igraph)
library(ggraph)
library(tidyverse)
library(scales)   # for pretty breaks

# 0. 读取边表（请确认文件路径/名称）
edges <- read.csv("ALS_Ligand_Receptor_Network.csv", stringsAsFactors = FALSE)

# 快速查看列名，确保列名是 ligand_ENSG / receptor_ENSG
print("Edges columns:")
print(colnames(edges))
# 示例期待: "ligand_ENSG" "receptor_ENSG" "log2FoldChange" "padj"

# 1. 规范列名（如果不是这个名字，重命名）
# 如果你的列确实是 ligand_ENSG / receptor_ENSG，这一步可以注释掉
if(!("ligand_ENSG" %in% colnames(edges)) & ("ligand_ensembl" %in% colnames(edges))){
  edges <- edges %>% rename(ligand_ENSG = ligand_ensembl,
                            receptor_ENSG = receptor_ensembl)
}

# 2. 清洗 ENSG 格式：去掉版本号、空格等（很关键）
clean_ensg <- function(x){
  x <- as.character(x)
  x <- ifelse(is.na(x), NA, x)
  x <- gsub("\\s+", "", x)        # 去空格
  x <- gsub("\\|.*", "", x)      # 去掉后缀如 "|GENE"
  x <- gsub("\\..*$", "", x)     # 去掉版本号如 .1
  return(x)
}

edges$ligand_ENSG  <- clean_ensg(edges$ligand_ENSG)
edges$receptor_ENSG<- clean_ensg(edges$receptor_ENSG)

# 3. 再次检查是否有 NA 或空行
cat("Edges rows:", nrow(edges), " Unique ligands:", length(unique(edges$ligand_ENSG)),
    " Unique receptors:", length(unique(edges$receptor_ENSG)), "\n")
cat("Any NA ligand:", any(is.na(edges$ligand_ENSG)), " Any NA receptor:", any(is.na(edges$receptor_ENSG)), "\n")

# 4. 为 igraph 构建边表（两个列名必须为 from/to）
# 建议先把列名改为 from / to 用于构图
edges_for_graph <- edges %>%
  select(from = ligand_ENSG, to = receptor_ENSG, log2FoldChange, padj)

# 5. 构建 igraph 对象（有向图）
g <- graph_from_data_frame(d = edges_for_graph, directed = TRUE)

# 6. 给节点打 type（Ligand / Receptor / Both）
V(g)$name[1:5]  # show examples

# 创建 sets
ligands_set   <- unique(edges_for_graph$from)
receptors_set <- unique(edges_for_graph$to)

V(g)$type <- "Other"
V(g)$type[ V(g)$name %in% ligands_set ]   <- "Ligand"
V(g)$type[ V(g)$name %in% receptors_set ] <- "Receptor"
# 若某节点同时是 ligand 和 receptor（存在双重身份），标为 "Both"
V(g)$type[ V(g)$name %in% intersect(ligands_set, receptors_set) ] <- "Both"

# 7. 映射 ligand 的 log2FC（note: receptors 可能没有值）
# 我们用 edges_for_graph 中的 from->log2FoldChange（若多个取平均）
ligand_fc_df <- edges_for_graph %>%
  group_by(from) %>%
  summarise(lfc = mean(log2FoldChange, na.rm = TRUE)) %>%
  as.data.frame()

V(g)$log2FC <- NA
match_idx <- match(V(g)$name, ligand_fc_df$from)
V(g)$log2FC[!is.na(match_idx)] <- ligand_fc_df$lfc[match_idx[!is.na(match_idx)]]

# 8. 诊断输出：看 type 分布和有无 log2FC
cat("Node counts by type:\n")
print(table(V(g)$type))
cat("Number of nodes with log2FC:", sum(!is.na(V(g)$log2FC)), "\n")

# 9. 计算 degree 并找 top hubs（degree 包括 in+out）
deg_all <- degree(g, mode = "all")
deg_sorted <- sort(deg_all, decreasing = TRUE)

topN <- 12  # 标注多少 hub
top_nodes <- names(head(deg_sorted, topN))
cat("Top hubs (by degree):\n")
print(head(deg_sorted, topN))

# 将 hub 名字记录到属性（方便绘图只标注 hub）
V(g)$label_hub <- ifelse(V(g)$name %in% top_nodes, V(g)$name, "")

# 10. 小网络提取（可选：只画 top ligands + 其 receptors，便于展示）
top_ligands <- names(sort(table(edges_for_graph$from), decreasing = TRUE))[1:20]
edges_small <- edges_for_graph %>% filter(from %in% top_ligands)
g_small <- graph_from_data_frame(edges_small, directed = TRUE)

# 11. 可视化（主图：全网络；小图：top20）
library(ggraph)
library(ggplot2)

# 全网络（只标注 top hubs）
p_all <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.4, edge_colour = "grey60") +
  geom_node_point(aes(color = type, size = pmax(0.5, abs(log2FC))), alpha = 0.9) +
  geom_node_text(aes(label = label_hub), repel = TRUE, size = 3, fontface = "bold") +
  scale_color_manual(values = c("Ligand" = "#E41A1C","Receptor" = "#377EB8","Both"="#4DAF4A","Other"="grey70")) +
  guides(size = guide_legend("abs(log2FC)")) +
  theme_void() +
  ggtitle("ALS Ligand–Receptor Network (Hub-labelled)")

# 小网络（更清晰，展示 top ligands）
p_small <- ggraph(g_small, layout = "fr") +
  geom_edge_link(alpha = 0.6, edge_colour = "grey40") +
  geom_node_point(aes(color = ifelse(name %in% edges_small$from, "Ligand","Receptor")),
                  size = 4) +
  geom_node_text(aes(label = ifelse(name %in% edges_small$from, name, "")), repel = TRUE, size = 3) +
  scale_color_manual(values = c("Ligand" = "#E41A1C","Receptor" = "#377EB8")) +
  theme_void() +
  ggtitle("Top-ligand centred LR subnetwork")

# 12. 显示图
print(p_all)
print(p_small)

# 13. 导出 network edge list + node attributes（方便 Cytoscape / DB）
# Edge list already exists; export node table
nodes_df <- data.frame(
  node = V(g)$name,
  type = V(g)$type,
  log2FC = V(g)$log2FC,
  degree = deg_all,
  label_hub = V(g)$label_hub,
  stringsAsFactors = FALSE
)
write.csv(nodes_df, "ALS_network_nodes.csv", row.names = FALSE)
write.csv(edges_for_graph, "ALS_network_edges_clean.csv", row.names = FALSE)

cat("Exported nodes and edges CSV. Done.\n")


