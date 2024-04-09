#!/usr/bin/env Rscript
# Be sure to first export RVERSION=4.3.2
# Also, have an ~/.Rprofile with .libPaths("/home/vhasfccuneod/Renv")

# install the R packages, if nbt already installed
for (package in c("plyr", "ggplot2", "segmented", "lcmm", "ggrepel", "parallel")) {
          if (!require(package, character.only=TRUE, quietly=TRUE)) {
                            install.packages(package, repos = "http://cran.us.r-project.org", lib="/home/vhasfccuneod/Renv")
    library(package, character.only=TRUE, quietly=TRUE)
          }
}
