# R scripts for analysis

To set up your R environment to run all of these scripts, run in the top
directory:

```
make boot
```

Note that all html files have `devtools::session_info()` called at the end,
which gives the R version along with the package versions. For packages
installed via GitHub, the specific commit address is documented and can be used
with the function `devtools::install_github()`.

## my\_functions.R

This R script contains various functions that are used with the following R
scripts. Documentation for each function is written in roxygen2's formatting
above each function. Functions intended for internal use have the keyword
"internal".

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
    isolates from Oregon and California](http://journals.plos.org/plospathogens/
    article?id=10.1371/journal.ppat.1000583).
 
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