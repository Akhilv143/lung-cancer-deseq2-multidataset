# Lung Cancer Bulk RNA-seq Meta-Analysis

Multi-cohort differential expression analysis of lung cancer integrating 3 GEO datasets.

## Datasets
| GEO ID | Platform | Samples |
|--------|----------|---------|
| GSE283245 | Illumina NovaSeq 6000 | Tumor + Control |
| GSE81089 | Illumina HiSeq 2500 | 19 matched pairs |
| GSE159857 | Illumina NextSeq 500 | LUAD + LUSC paired |

## Workflow
1. Metadata harmonization & count matrix merging
2. ComBat-seq batch correction
3. DESeq2 differential expression (Tumor vs Control)
4. GO (BP/CC/MF) and KEGG enrichment analysis
5. Heatmap, Volcano plot, PCA visualization

## Key Results
- Significant DEGs (padj < 0.05, |log2FC| > 1): see `results_LUNG/02_tables/DEG/`
- Enrichment results: see `results_LUNG/02_tables/enrichment/`

## Requirements
R packages: DESeq2, sva, clusterProfiler, ComplexHeatmap, ggrepel, org.Hs.eg.db
