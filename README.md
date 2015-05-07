# Supplementary Material for poppr 2.0 manuscript

This repository contains supplemental materials relating
to the manuscript submitted to Frontiers called "Novel R tools for analysis of genome-wide population genetic data with emphasis on clonality" By Zhian Kamvar (@zkamvar), Jonah Brooks (@JonahBrooks), and Nik Grunwald (@grunwald).

# Citation

Cite this supplementary material as:

> ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementray Material for Frontiers Plant Genetics and
Genomics 'Novel R tools for analysis of genome-wide population genetic data with emphasis on clonality'. DOI: XXX.

# Setup

To use the scripts in this repository, you will need the latest versions of poppr and adegenet. 

### From a unix shell:

```sh
make boot
```
### From within R:

```R
local({r <- getOption("repos"); r["CRAN"] <- "http://cran.at.r-project.org"; options(repos = r)})
install.packages(c("ape", "knitr", "rmarkdown", "animation", "devtools"))
devtools::install_github(c("thibautjombart/adegenet", "emmanuelparadis/pegas/pegas", "KlausVigo/phangorn"))

devtools::install_github("grunwaldlab/poppr@2.0-rc")
```

# Folders

## Rscripts

The folder called **Rscripts** contain R scripts, corresponding HTML files, and cached information used for the manuscript and supplementary information. Please see `Rscritps/README` for more details.

## images

The folder called **images** contains HTML animations of dendrograms collapsing multilocus genotypes by genetic distance. 
