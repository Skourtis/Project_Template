#Settting up directory
#install.packages("tidyverse")
pacman::p_load(piggyback, renv, here, tidyverse )
#renv::restore()

## Created a first release on Github
#pb_new_release("Skourtis/Project_Template")
piggyback::pb_track("*.zip") %>%
    pb_upload(repo = "Skourtis/Project_Template")

##end
renv::snapshot()