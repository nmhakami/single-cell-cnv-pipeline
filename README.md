# single-cell-cnv-pipeline
A pipeline for generating and analyzing single cell copy number profiles from whole genome sequencing data.

## Overview
This pipeline processes sequencing data to generate genome-wide copy number profiles and perform downstream quality control and clustering. 
The workflow consists of three main stages: 
1) Sequencing read processing and copy number profile generation
2) Quality assessment and filtering of noisy CNV profiles
3) HDBSCAN-based clustering of high-quality CNV profiles to categorize cells as flat, clonal, or rearranged but not clonal.

## Pipeline
FASTQ -> Alignment -> Read Filtering -> Bin Counting -> GC Correction -> Segmentation -> CNV Profile -> QA/QC -> HDBSCAN clustering
   
