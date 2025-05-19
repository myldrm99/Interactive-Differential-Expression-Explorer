# Interactive Differential Expression Explorer

[![Shiny App Screenshot](Screenshot%202025-05-17%20at%2020.16.08.png)](https://muhammedyildidirm.shinyapps.io/mbg513/)

A Shiny web application for interactive analysis of GEO datasets and user-uploaded expression data with advanced differential-expression visualization.

---

## 🎥 Tutorial Video

Watch the five-minute walkthrough:  
[![Play tutorial_4.mp4](https://img.shields.io/badge/Tutorial-▶️%20Play-green.svg)](tutorial_4.mp4)  
*(Click above or navigate to `tutorial_4.mp4` in this repo.)*

---

## 🚀 Key Features

- **Real-time analysis** with adjustable thresholds (log₂ FC & adjusted *p*-value)  
- **Interactive sample grouping** via clickable data table  
- **Rich visualizations**: volcano plots, MA/scatter plots, heatmaps, PCA  
- **Dual data sources**: GEO accession IDs or custom CSV/TSV uploads  
- **One-click exports**: download all figures (PNG) & results tables (CSV)

---

## 🛠 Installation

```r
# Install required packages from CRAN
install.packages(c(
  "shiny",
  "GEOquery",
  "DT",
  "limma",
  "ggplot2",
  "plotly",
  "pheatmap"
))
