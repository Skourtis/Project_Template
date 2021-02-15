#Settting up directory
##Settting up directory
install.packages("pacman")
pacman::p_load(piggyback, renv, here, tidyverse )


## Created a first release directly on Github
#pb_new_release("Skourtis/Project_Template")
piggyback::pb_track("*.zip") %>%
    pb_upload(repo = "Skourtis/Project_Template")

##end
renv::snapshot()