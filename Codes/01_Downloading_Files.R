#Downloading from Git template
options(repos = getOption("repos")["CRAN"])
renv::restore()
devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)
pacman::p_load(piggyback, renv, here, tidyverse)

pb_download()
