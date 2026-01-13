# Spatial Analysis of Tumor–Macrophage Niches Using Visium HD Data

This repository contains a fully reproducible spatial transcriptomics analysis pipeline applied to **raw Visium HD (8 µm) data** from a human colon cancer sample.  
The goal of the project is to independently reproduce and validate key spatial findings reported in a published study by directly analyzing raw spatial data and generating quantitative, statistical, and visual outputs.

---

## Project Overview

Spatial transcriptomics enables the investigation of how gene expression programs are organized in physical space. In colorectal cancer, the tumor boundary (tumor periphery) represents a biologically active interface where immune and stromal interactions are spatially structured.

In this project, we:

- Define tumor regions using **CEACAM5/CEACAM6** expression
- Construct a **distance-based tumor periphery** (≤ 50 µm)
- Identify **macrophage-enriched spatial bins** using **C1QC**
- Classify periphery macrophages into **SPP1⁺** and **SELENOP⁺** subtypes
- Perform **non-parametric statistical comparisons**
- Generate **publication-ready figures and supplementary tables**

All analyses are performed directly from **raw Visium HD outputs**, without relying on cell segmentation or pre-annotated regions.

---

## Dataset

- **Technology:** 10x Genomics Visium HD  
- **Resolution:** 8 µm binned spatial transcriptomics  
- **Sample:** Human colon cancer (P1)  
- **Input files include:**
  - `filtered_feature_bc_matrix.h5`
  - `tissue_positions.parquet`
  - Full-resolution spatial coordinates

The raw dataset is publicly available from 10x Genomics.

---

## Analysis Pipeline

The spatial analysis workflow consists of the following steps:

1. **Data loading**
   - Read expression matrix from `filtered_feature_bc_matrix.h5`
   - Load spatial coordinates from `tissue_positions.parquet`

2. **Normalization**
   - Library-size normalization (counts per 10,000)
   - Log1p transformation

3. **Tumor region identification**
   - Tumor score = CEACAM5 + CEACAM6
   - Tumor bins defined as top 10% of the score distribution

4. **Tumor periphery definition**
   - Euclidean distance to nearest tumor bin
   - Periphery defined as bins within 50 µm

5. **Macrophage identification**
   - Macrophage-enriched bins defined by C1QC expression

6. **Macrophage subtype classification**
   - SPP1⁺ vs SELENOP⁺ based on relative marker expression

7. **Statistical analysis**
   - Mann–Whitney U test for subtype comparisons
   - Median expression and p-values reported

8. **Visualization and reporting**
   - Spatial maps
   - Distance distributions
   - Gene expression boxplots
   - Excel-based supplementary tables

---

## Repository Structure

```text
.
├── spatial_report.py          # Main analysis script
├── requirements.txt           # Python dependencies
├── P1_8um_min.tar.gz           # Raw Visium HD input data
├── outputs/
│   ├── P1_spatial_report_outputs.xlsx
│   ├── Fig1_tumor_periphery.png
│   ├── Fig2_macrophage_subtypes.png
│   ├── Fig3_distance_hist.png
│   ├── Fig4_key_gene_boxplots.png
│   └── run_summary.json
└── README.md
```

## Outputs

Excel (Supplementary Table)

P1_spatial_report_outputs.xlsx contains:

- Region-level bin counts (tumor, periphery, tissue)

- Macrophage prevalence by region

- Macrophage subtype counts

- Median distance-to-tumor statistics

- Mann–Whitney U test results

- Full bin-level spatial annotation table

This file enables full reproducibility and independent re-analysis.

## Figures

- Figure 1: Tumor and tumor periphery spatial map

- Figure 2: Spatial distribution of SPP1⁺ and SELENOP⁺ macrophages

- Figure 3: Distance-to-tumor distributions

- Figure 4: Differential gene expression between macrophage subtypes

All figures are generated directly from raw spatial data.

## Key Results (Summary)

- Total spatial bins analyzed: 507,684

- Tumor bins: 50,721

- Tumor periphery bins: 123,874

- Periphery macrophage bins: 2,769

- SPP1⁺ macrophages: 54

- SELENOP⁺ macrophages: 2,715

These results demonstrate strong spatial structuring of macrophage programs at the tumor boundary.

## Reproducibility

- No manual region annotation was used

- All thresholds are explicitly defined

- All results are traceable to Supplementary Tables

- The pipeline can be applied to other Visium HD datasets with minimal modification

## Scientific Context

This project was developed as part of a graduate-level spatial transcriptomics assignment and serves as an independent validation of published findings using raw Visium HD data.

It demonstrates how biologically meaningful spatial patterns can be recovered from bin-level data without requiring high-level proprietary tools.