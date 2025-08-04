# EdgeR

This repository contains a Jupyter notebook demonstrating an RNA-seq analysis pipeline using the **edgeR** and **limma** Bioconductor packages. The notebook installs required dependencies, imports transcript quantification results from *kallisto*, performs quality control, normalizes counts with TMM, and calculates differential expression statistics with the limma-voom workflow.

The analysis is run from Python via the `rpy2` extension so that R code can be executed inside notebook cells. Sample metadata are joined to the count matrix and a design matrix is generated to test group-wise contrasts.

The notebook also includes utility plots such as log-CPM boxplots and IQR versus median checks to assess sample quality.

## Colab interactive usage

Run the pipeline interactively in Google Colab by executing `%run EdgeR_colab.py` in a notebook cell. Use the displayed widgets to
select the metadata columns corresponding to sample names and groups, then click **Run EdgeR pipeline** to launch the analysis
with those selections.
