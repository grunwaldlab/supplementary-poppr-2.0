# Supplementary Material for poppr 2.0 manuscript

This repository contains supplemental materials relating
to the manuscript submitted to Frontiers called "Novel R tools for analysis of genome-wide population genetic data with emphasis on clonality" By Zhian N. Kamvar ([@zkamvar](https://github.com/zkamvar)), Jonah C. Brooks ([@JonahBrooks](https://github.com/JonahBrooks)), and Niklaus J. Grünwald ([@grunwald](https://github.com/grunwald)).

# Citation

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.17424.svg)](http://dx.doi.org/10.5281/zenodo.17424)

Please cite this supplementary material as:

> ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary Material for Frontiers Plant Genetics and
Genomics 'Novel R tools for analysis of genome-wide population genetic data with emphasis on clonality'.  Zenodo. [10.5281/zenodo.17424](http://dx.doi.org/10.5281/zenodo.17424)

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

## Images

The folder called **images** contains HTML animations of dendrograms collapsing multilocus genotypes by genetic distance. 
