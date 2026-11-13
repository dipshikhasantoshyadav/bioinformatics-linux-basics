
# SARS\\oV-2 RNA-seq Differential Gene Expression Analysis

## Project Overview
This project investigates transcriptional changes in normal human bronchial epithelial (NHBE) cells.

The objective is to identify differentially expressed genes between mock-treated control and SARS-CoV-2-infected samples and explore the biological processes associated with the observed changes.

## Dataset
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Accession:** GSE147507
- **Samples:** 6 NHBE samples
- **Groups:** 3 mock controls and 3 SARS-CoV-2-infected samples
- **Comparison:** Infected vs Control

## Analysis Workflow
1. Retrieved RNA-seq read count data from GEO.
2. Selected the six NHBE samples.
3. Prepared the gene count matrix and sample metadata.
4. Performed differential gene expression analysis using DESeq2.
5. Identified significantly upregulated and downregulated genes.
6. Performed Gene Ontology Biological Process enrichment using clusterProfiler.
7. Visualised the results using GO enrichment and volcano plots.

## Key Findings
- 91 significantly upregulated genes were identified using adjusted p-value < 0.05 and absolute log2 fold change > 1.
- GO enrichment highlighted inflammatory responses, cytokine production, viral defence and chemotaxis.
- Notable inflammatory genes included IL6, IL1A, IL1B, TNF, CXCL2, CXCL5 and CXCL6.

- 14 significantly downregulated genes were identified using adjusted p-value < 0.05 and absolute log2 fold change > 1.


## Tools and Packages
- R
- DESeq2
- clusterProfiler
- org.Hs.eg.db
- ggplot2

## Repository Structure
- `data/` – GEO count data, sample metadata and analysis outputs
- `scripts/` – Analysis scripts
- `results/` – Final figures and results
- `documentation/` – Project notes and methodology

## Reproducibility
The analysis uses publicly available GEO data. Sample metadata and processed count data are included to document the workflow.

## Author
Dipshikha Yadav

