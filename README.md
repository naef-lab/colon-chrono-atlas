# Chrono-atlas of cell-type specific daily gene expression rhythms in the undamaged and regenerating mouse colon

## Overview

This repository contains the analysis code for the single-cell RNA-seq chrono-atlas of the mouse colon, covering undamaged and DSS-induced regenerating conditions. The study profiled 24-hour gene expression rhythms across >20 cell types in the proximal and distal colon.

## Data availability

| Resource | Accession / DOI | Description |
|----------|----------------|-------------|
| ENA | [PRJEB102541](https://www.ebi.ac.uk/ena/browser/view/PRJEB102541) | Raw sequencing reads (FASTQ) |
| Zenodo | [DOI: XXXX](https://doi.org/XXXX) | Cell Ranger raw h5, CellBender-corrected h5, final annotated h5ad |

## Pipeline overview

The analysis proceeds in six steps:

1. **Cell Ranger** (v7.0.1): alignment to GRCm39 and quantification (10x Flex assay)
2. **CellBender** (remove-background): ambient RNA removal, expected cells = 5,000
3. **DoubletFinder**: doublet detection at 4% rate with pK optimization
4. **Preprocessing** (Seurat v4.1.3): QC filtering, rPCA integration, Leiden clustering, cell type annotation, ISC state scoring
5. **Rhythmicity analysis**: pseudo-bulk harmonic regression (24h period), cSVD decomposition
6. **Pathway enrichment**: enrichR (KEGG, Reactome, WikiPathways)

## Repository structure

```
colon-chrono-atlas/
├── 01_cellbender/          # Ambient RNA removal (SLURM + CellBender)
├── 02_doubletfinder/       # Doublet detection (SLURM + R)
├── 03_preprocessing/       # QC, integration, clustering, annotation
├── 04_rhythmicity/         # Harmonic regression and cSVD
├── 05_figures/             # Figure generation code
└── 06_enrichment/          # Pathway enrichment analysis
```
