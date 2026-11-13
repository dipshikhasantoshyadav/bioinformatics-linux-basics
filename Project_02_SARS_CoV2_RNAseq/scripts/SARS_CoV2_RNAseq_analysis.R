# =====================================================
# SARS-CoV-2 RNA-seq Differential Expression Analysis
# Dataset: GSE147507
# Comparison: Infected vs Control (NHBE cells)
# Author: Dipshikha Yadav
# =====================================================

# 1. Load required packages

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 2. Load gene count data

counts_data <- read.delim(
  "data/NHBE_SARSCoV2_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

# 3. Load sample metadata

sample_info <- read.csv(
  "data/sample_metadata.csv",
  row.names = 1
)

# 4. Set the condition as a factor

sample_info$condition <- factor(
  sample_info$condition,
  levels = c("Control", "Infected")
)

# 5. Check the data

dim(counts_data)
head(counts_data)
sample_info

# 6. Create DESeq2 dataset

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_data)),
  colData = sample_info,
  design = ~ condition
)

# 7. Set Control as reference

dds$condition <- relevel(dds$condition, ref = "Control")

# 8. Run differential expression analysis

dds <- DESeq(dds)

# 9. Extract results

res <- results(dds)

# 10. View top results

head(res)

# 11. Filter significant genes

sig <- res[
  !is.na(res$padj) &
  res$padj < 0.05 &
  abs(res$log2FoldChange) > 1,
]

# 12. Identify upregulated genes

upregulated <- sig[sig$log2FoldChange > 1, ]

# 13. Identify downregulated genes

downregulated <- sig[sig$log2FoldChange < -1, ]

# 14. Count significant genes

nrow(upregulated)
nrow(downregulated)

# 15. Save differential expression results

write.csv(as.data.frame(res),
          "data/DESeq2_all_results.csv")

write.csv(as.data.frame(upregulated),
          "data/Upregulated_genes.csv")

write.csv(as.data.frame(downregulated),
          "data/Downregulated_genes.csv")

# 16. Convert gene symbols into Entrez IDs

gene_ids <- bitr(
  rownames(upregulated),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 17. Perform GO Biological Process enrichment

ego_up <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 18. Convert GO results into a table

go_table <- as.data.frame(ego_up)

# 19. Save GO results

write.csv(
  go_table,
  "data/GO_upregulated_BP_results.csv",
  row.names = FALSE
)


# 20. Generate GO enrichment dot plot

go_plot <- dotplot(ego_up, showCategory = 20)

# 21. Save GO enrichment plot

ggsave(
  "data/GO_enrichment_dotplot.png",
  plot = go_plot,
  width = 10,
  height = 8,
  dpi = 300
)


# 22. Prepare data for volcano plot

res_df <- as.data.frame(res)

res_df$group <- "Not significant"

res_df$group[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange > 1
] <- "Upregulated"

res_df$group[
  !is.na(res_df$padj) &
  res_df$padj < 0.05 &
  res_df$log2FoldChange < -1
] <- "Downregulated"


# 23. Generate volcano plot

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = group
  )
) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "SARS-CoV-2 vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Category"
  )


# 24. Save volcano plot

ggsave(
  "data/SARS_CoV2_volcano_plot_coloured.png",
  plot = volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)


