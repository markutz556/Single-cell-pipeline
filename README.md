# Analysis Pipeline for Fluidigm Polaris Platform

<br></br>
## Table of Contents
1. [Fluidigm Polaris](#fluidigm-polaris)
2. [Data source](#data-source)
3. [Code structure](#code-structure)
4. [Publication](#publication)
5. [Acknowledgement](#acknowledgement)

<br></br>
## Fluidigm Polaris platform
 * [Fluidigm Polaris platform](https://www.fluidigm.com/products/polaris) is a full length single-cell analysis technology that allows to actively select single cells from heterogeneous cell samples and hold them in individual, environmentally controlled reaction chambers. It can correlate a cell’s phenotype, physiological condition and molecular signature to discover the genetic circuits regulating cell function or compare a population’s sensitivity to environmental changes. 

<br></br>
## Data source
 * Erythroleukemia [K562](https://www.ncbi.nlm.nih.gov/pubmed/163658) cell line as control. 
 * Breast cancer [SUM149](https://www.ncbi.nlm.nih.gov/pubmed/17157791) cell line. 
 * Preprocessing of single cells in [Fluidigm Polaris platform](https://www.fluidigm.com/products/polaris).
 * Sequencing on [Illumina HiSeq 2500](https://www.illumina.com/systems/sequencing-platforms/hiseq-2500.html).
 * Analysis using [R](https://www.r-project.org).
   
<br></br>
## Code structure
1. Import data
2. Density plot
3. Differential expression
4. Volcano plot
5. GO analysis
6. PCA
7. Visualization (bar plot, box plot, dot plot, violin plot)

<br></br>
## Publication
 <I>Analytic Pipeline for full-length single-cell RNA sequencing data acquired by Fluidigm Polaris sequencing platform</I>

<br></br>
## Acknowledgement
 * Guan lab, University of Michigan, Ann Arbor
 * DCMB, University of Michigan, Ann Arbor
