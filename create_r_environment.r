#!/usr/bin/env Rscript
# Environment variables must be setup first.
##
# export RVERSION=4.3.2
# export R_LIBS_USER=~/Renv


# The library `parallel`, seems to be included in the base distribution and throws a warning when trying to install.
# install the R packages, if nbt already installed
for (package in c("plyr", "ggplot2", "segmented", "lcmm", "ggrepel", "reshape2", "devtools")) {
          if (!require(package, character.only=TRUE, quietly=TRUE)) {
                            install.packages(package, repos = "http://cran.us.r-project.org", dependencies=TRUE, lib="/home/vhasfccuneod/Renv")
    library(package, character.only=TRUE, quietly=TRUE)
          }
}
