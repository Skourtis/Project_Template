#Settting up directory
install.packages("pacman")
pacman::p_load(piggyback, renv, here, tidyverse)

piggyback::pb_track(here::here("Datasets")) %>%
    pb_upload()

##end
renv::snapshot()
