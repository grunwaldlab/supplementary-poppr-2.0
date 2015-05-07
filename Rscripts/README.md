# R scripts for analysis

To set up your R environment to run all of these scripts, run in the top
directory:

```
make boot
```

Note that all html files have `devtools::session_info()` called at the end,
which gives the R version along with the package versions. For packages
installed via GitHub, the specific commit address is documented and can be used
with the function `devtools::install_github()`. In the manuscript, the *poppr*
version is stated as being version `2.0` whereas at the bottom of all of these
html files, the version will appear as `1.1.4.99.x`. This is intentional. 
`1.1.4.99` is to signify that this is a development version after version `1.1.4`
and the number that comes after 99 indicates the number of changes that have
occurred since verion `1.1.4`. The version will increase to 2.0 once it has
been submitted to CRAN.

## my\_functions.R

This R script contains various functions that are used with the following R
scripts. Documentation for each function is written in roxygen2's formatting
above each function. Functions intended for internal use have the keyword
"internal".

To use these functions without downloading the repository, you can use devtools:

```R
devtools::source_url("https://raw.githubusercontent.com/grunwaldlab/supplementary-poppr-2.0/master/Rscripts/my_functions.R")
```

## figure\_one.R

Construction of figure 1 with a walkthrough of the concepts using example data.

## mlg\_filter\_stats.R

This script's purpose is to explore how the MLG filtering methods perform in two
ways:

 1. With *P. infestans* data from Genotype-ID
 2. With simulated data

## msn\_graph\_walk.R

The purpose of this script is to demonstrate the differences between minimum
spanning trees with and without reticulation. Two data sets are used:

 1. *Phytophthora ramorum* from the [Oregon Sudden Oak Death
    Epidemic](http://dx.doi.org/10.5281/zenodo.13007) combined with [Nursery
    isolates from Oregon and California](http://www.plospathogens.org/article/info%3Adoi%2F10.1371%2Fjournal.ppat.1000583).
 
 2. Sequences from the HA segment of the H3N2 flu virus from the package
    [adegenet](https://github.com/thibautjombart/adegenet).
 
Both of these data sets produce minimum spanning networks with very few
branching events without reticulation. With reticulation, they create dense
clusters.

We use the [infoMAP](http://www.pnas.org/content/105/4/1118)
[algorithm](http://arxiv.org/abs/0707.0609) to detect clusters within the graphs
and then compared the entropy.

## index\_of\_association.R

This file simulates data sets that contain population structure and also
simulates clonal data sets. These are compared with panmictic data sets with a
sliding window of the index of association and random sampling.
