# guenon_phylogenomics
Custom python scripts used in Jensen et al. 2023.

## Listing the scripts used per subsection:
### Whole genome phylogenetic analyses:
- windowPhylogenies.py
  Used to build maximum likelihood and NJ trees in sliding windows (Using IQTree and phyML). Modified from Simon Martins phyml_sliding_window.py (https://github.com/simonhmartin/genomics_general)
- treeParser.py
  Various parsing and pruning of trees.

### Gene flow analyses
- findFixedDifferences.py
  Used to find SNPs differentially fixed between two specified groups of samples. Used in the private D-statistics and N-mt analyses.
- classifyTopology.py
  Used to classify window topologies into one of max 5 categories as described in the manuscript.

### Mitochondrial phylogeny
The following scripts were used to parse and process the mitochondrial assemblies, prior to phylogenetic analyses:
