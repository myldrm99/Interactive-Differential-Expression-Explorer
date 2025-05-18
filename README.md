# Interactive Differential Expression Explorer

[![Shiny App Screenshot](Screenshot%202025-05-17%20at%2020.16.08.png)](https://muhammedyildidirm.shinyapps.io/mbg513/)

A Shiny web application for interactive analysis of GEO datasets and user-uploaded expression data with advanced differential expression visualization.

## Key Features
- **Real-time analysis** with adjustable thresholds (log FC & p-value)
- **Interactive sample grouping** via clickable data table
- **Multiple visualizations**: Volcano plots, MA plots, Heatmaps, PCA
- **Dual data sources**: GEO accession IDs or custom CSV/TSV uploads
- **One-click exports** for all figures and results tables

## Installation
```r
# Required R packages
install.packages(c("shiny", "GEOquery", "DT", "limma", "ggplot2", "plotly", "pheatmap"))
```

## Usage
1. Run app locally:
```r
shiny::runGitHub("yourusername/your-repo-name")
```
2. Choose data source (GEO ID or file upload)
3. Select samples using table checkboxes
4. Adjust analysis thresholds using sliders
5. Explore results in interactive tabs
6. Download publication-ready figures/results

## Dependencies
- Core: `shiny`, `GEOquery`, `limma`
- Visualization: `ggplot2`, `plotly`, `pheatmap`
- UI: `DT`

*Report and code available in repository*  
📈 [Live Demo](https://muhammedyildidirm.shinyapps.io/mbg513/) | 📚 [GEO Database](https://www.ncbi.nlm.nih.gov/geo/)
