# sensoryfates_paper

Code for reproducibiliy

required forked packages:
- https://github.com/LouisFaure/crestree
- https://github.com/LouisFaure/pagoda2 

## Logic of the analysis

1. [QC and preprocessing](Preprocessing_notebook.md) of raw count matrices in R using pagoda2 R package. 
2. [RNA velocity](Velocity_notebook.ipynb) analysis of raw spliced/unspliced count matrices in python using scvelo package. Mapping velocity on UMAP generated during the preprocessing part.
3. [Pseudotime](Pseudotime_analysis.md) analysis, including count matrix correction with scde R package, tree fitting using ElPiGraph.R and downstream pseudotime analysis using crestree R package

## Citation

Faure, L., Wang, Y., Kastriti, M.E. et al. Single cell RNA sequencing identifies early diversity of sensory neurons forming via bi-potential intermediates. Nat Commun 11, 4175 (2020). https://doi.org/10.1038/s41467-020-17929-4
