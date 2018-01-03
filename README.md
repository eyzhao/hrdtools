# HRDtools

HRDtools is an R package which computes HRD scores and mutation signatures in individual tumours.

## Installation

HRDtools is an R package, and can be installed via devtools.

```{r}
devtools::install_github('eyzhao/hrdtools')
```

HRDtools has numerous package dependencies, which may need to be installed first.

```{r}
install.packages(c(
  'tidyverse',
  'nnls',
  'parallel',
  'snow',
  'snowfall',
  'plyr'
))

source("https://bioconductor.org/biocLite.R")
bioclite(c(
  'GenomicRanges',
  'VariantAnnotation',
  'SomaticSignatures',
  'BSgenome'
))
```

## Usage

Once installed, HRDtools can be used in R as follows,

```{r}
library(hrdtools)

run_test(
  'path/to/segments.file.tsv'
)
```

where the file `path/to/segments.file.tsv` contains CNV/LOH segments. HRDtools was tested on output from [APOLLOH](http://shahlab.ca/projects/apolloh/), but in theory any CNV caller with LOH support should work.

The segments file should contain one segment per row, with the columns `chr, start, end, copynumber, lohtype`. If different column names are used to denote copy number and LOH, they can be specified as follows.

```{r}
run_test(
  'path/to/segments.file.tsv',
  loh.col = 'loh_state',
  cnv.col = 'tumour_copy_number'
)
```

## Citation

If you use HRDtools in your publication, please cite the following study:

[Zhao EY, Shen Y, Pleasance E, ... Jones SJM, Homologous Recombination Deficiency and Platinum-Based Therapy Outcomes in Advanced Breast Cancer. Clinical Cancer Research. 2017 Dec 15;23(24):7521-7530](https://www.ncbi.nlm.nih.gov/pubmed/29246904)
