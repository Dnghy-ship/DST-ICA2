
library(tidyverse)
library(biomaRt)

#1. 读入 Ramilowski 配体–受体数据库

lr_raw <- read.delim("PairsLigRec.txt", header = TRUE, sep = "\t")

# 先统一列名
colnames(lr_raw)[2] <- "ligand_symbol" 
colnames(lr_raw)[4] <-"receptor_symbol"

# 2. 分别提取所有 ligand 和 receptor symbol

ligands <- unique(lr_raw$ligand_symbol)
receptors <- unique(lr_raw$receptor_symbol)

all_symbols <- unique(c(ligands, receptors))

# 3. 用 biomaRt 将 symbol 转成 ENSG
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = all_symbols,
  mart = mart
)
# 去掉空行
mapping <- mapping %>% filter(ensembl_gene_id != "")

# 4. 把 mapping merge 回 L–R 表
# merge ligand
lr_mapped <- lr_raw %>%
  left_join(mapping, by = c("ligand_symbol" = "hgnc_symbol")) %>%
  rename(ligand_ENSG = ensembl_gene_id) %>%
  left_join(mapping, by = c("receptor_symbol" = "hgnc_symbol")) %>%
  rename(receptor_ENSG = ensembl_gene_id)
# 去掉没有 ENSG 的行
lr_mapped <- lr_mapped %>%
  filter(!is.na(ligand_ENSG), !is.na(receptor_ENSG))

# 5. 读入 secreted proteins ENSG 列表
secreted <- read.table("Secreted_protein_ref_ENSEMBLE_union.txt", stringsAsFactors = FALSE)[,1]

# 6. 过滤：只保留 secreted 作为 ligand 的 pairs
result <- lr_mapped %>%
  filter(ligand_ENSG %in% secreted) %>%
  select(ligand_symbol, ligand_ENSG,
         receptor_symbol, receptor_ENSG) %>%
  distinct()

# 7. 输出结果
write.csv(result, "Secreted_to_Receptors_Ramilowski.csv", row.names = FALSE)

print("处理完成！输出文件：Secreted_to_Receptors_Ramilowski.csv")
