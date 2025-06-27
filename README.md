# 📂 R

This folder contains custom R packages developed to analyze various types of omics data. Each package is designed to be modular and focused on specific analytical tasks, enabling streamlined and reproducible workflows for bioinformatics projects.

---

## 📦 Contents

### ✅ `genomicsFun`

**Description:**  
A set of R functions for the quality control and visualization of mutational data from whole-exome sequencing (WES) and whole-genome sequencing (WGS).  

**Main features:**  
- Sample-level QC metrics
- Mutation spectrum plots
- Mutation burden estimation
- Context-specific mutational signatures
- Publication-ready visualizations

**Intended use:**  
This package supports researchers working with cancer genomics, germline studies, or any project that requires exploratory and QC analysis of variant call data from WES/WGS pipelines.

---

## 💡 General structure

Each package in this folder is designed to be:

- **Self-contained**: All dependencies and functions are documented.
- **Reproducible**: Includes example data or minimal reproducible workflows.
- **Extensible**: Easy to adapt to new data types (e.g. methylation, transcriptomics, proteomics).

---

## ⚙️ Future packages (planned)

- `methylationFun`: Tools for processing and visualizing DNA methylation array or WGBS data.
- `transcriptomicsFun`: QC and differential expression analysis for RNA-seq data.

---

## 📜 How to install

You can install these packages (example below for `genomicsFun`) using:

```r
# Example
devtools::install_local("R/genomicsFun")
```

## 📜 How to load
```r
# Example
devtools::load_all("R/genomicsFun")
```