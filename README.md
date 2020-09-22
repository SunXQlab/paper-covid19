# Tissue Microenvironment-mediated Expression of SARS-CoV-2 Receptor ACE2

We collected scRNA-Seq data and bulk RNA-Seq data from COVID-19 patients BALF to analyze microenvironment-mediated ACE2 regulation after SARS-CoV-2 infection. We developed a single-cell transcriptomic data-based multi-layer signaling network approach to systematically construct the inter-/intracellular signaling pathways regulating ACE2 expression based on the integration of prior information, gene expression and statistical tests.   

## Workflow

The codes were implemented in R (version 3.6.1). The codes included two parts:

1. covid_balf_sc (scRNA-Seq data)
    The `/covid_balf_sc/`  folder includes codes and data in analysis of scRNA-Seq data from COVID-19 patients BALF, codes and databases of scMLnet algorithm and genesets of target gene in comparison between scMLnet and NicheNet.

2. covid_balf_bulk (bulk RNA-Seq data)
    The `/covid_balf_bulk/`  folder includes codes, data and databases in analysis of bulk RNA-Seq data from COVID-19 patients BALF.
