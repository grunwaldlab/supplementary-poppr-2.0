boot: inst_deps ver_two_deps inst_poppr_two 

inst_deps:
	R --slave -e 'install.packages(c("ape", "knitr", "rmarkdown", "animation", "devtools"), repos = "http://cran.at.r-project.org")'; \

inst_poppr_one:
	R --slave -e "local({r <- getOption("repos"); r["CRAN"] <- "http://cran.at.r-project.org"; options(repos = r)}); devtools::install_github('grunwaldlab/poppr')"; \

inst_poppr_two:
	R --slave -e "local({r <- getOption("repos"); r["CRAN"] <- "http://cran.at.r-project.org"; options(repos = r)}); devtools::install_github('grunwaldlab/poppr@2.0-rc')";

ver_two_deps:
	R --slave -e 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.at.r-project.org"; options(repos = r)}); devtools::install_github(c("thibautjombart/adegenet", "emmanuelparadis/pegas/pegas", "KlausVigo/phangorn"))';

