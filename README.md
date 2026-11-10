# Spatial transcriptomics reveals inflammatory fibro-adipogenic progenitors within fibro-inflammatory niches in Duchenne muscular dystrophy 
This repository contains the analysis scripts for the results in the above paper.

# Data
The processed files for snRNA-seq reference are available at- https://doi.org/10.5281/zenodo.20718391
The raw and processed files for the Visium dataset are available at GEO under the accession number GSE353421 (publicly available from 16th Nov 2026).

# Requirements

## Python (3.10.18)
- anndata 0.11.4
- cell2location 0.1.5
- liana 1.7.1
- matplotlib 3.10.3
- numpy 1.26.4
- pandas 2.3.0
- scanpy 1.11.2
- seaborn 0.13.2
- squidpy 1.6.5

## R (4.5.3)
- Seurat 5.3.0
- DoubletFinder 2.0.6
- Harmony 1.2.3

# Scripts
1.snRNA_faps_subtyping.R - Seurat based workflow for generating DMD with FAPs subtypes final object that is used as reference for the deconvolution.
2.Visium_preprocessing.py- Pre-processing of 10X visium raw data of DMD samples.
3.Visium_cell2location_deconvolution.py - Cell2location deconvolution workflow
4.Visium_NMF_neighborhoodenrichment - NMF analysis and neighbourhood enrichment of the niches.
5.Visium_CC_communication1.py - LIANA+ Misty workflow for communication analysis between niches and cell types.
6.Visium_CC_communication2.R - Consolidating resulting generated in Visium_CC_communication1.py
7.CC_visualisation.R - Generating plots for visualisation of the CC communication analysis results.
8.Niche_functionalenrichment.R - Pathway enrichment of the DEGs identified for the niches.


