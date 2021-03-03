### following DEP
BiocManager::install("DEP")
install.packages("janitor")
library("janitor")
read_tsv(paste0(txt_folder,
       "/proteinGroups.txt"))
data <- UbiLength
protein_groups <- read_tsv(paste0(txt_folder,
                                  "/proteinGroups.txt")) %>% janitor::clean_names()
data <- subset(protein_groups, is.na(protein_groups$reverse == "+") &
                   is.na(protein_groups$potential_contaminant == "+")) 
experimental_design  <- data.frame(label = 1:9,
                      condition = c(rep("Control",3),
                                    rep("sh1",3),
                                    rep("sh2",3)),
                      replicate = rep(1:3,3))
LFQ_columns <- grep("lfq_", colnames(data))
data_unique <- make_unique(data, "majority_protein_i_ds", "majority_protein_i_ds", delim = ";")
data_se <- make_se(data_unique, LFQ_columns, experimental_design)




n_of_replicates <- 3
condition <- c("Control","sh1","sh2")
data <- data %>% dplyr::select(contains("majority")|contains("lfq_")) %>%
    set_names(c("Uniprot",paste0("Control_",1:n_of_replicates), paste0("sh1_",1:n_of_replicates),paste0("sh2_",1:n_of_replicates))) %>%
    column_to_rownames("Uniprot") %>% 
    mutate(across(.cols = everything(),log10))
    
data_per_condition <- purrr::map(.x = condition,~dplyr::select(data,contains(.x)))

imputation_proteomics <- function(x){
    x = data[[1]]
    deciles <- quantile(x %>% unlist() %>% subset(.,.>0), #### finding the bottom 2 deciles of the distribution 
             prob = seq(0, 1, length = 11), type = 5) %>%
        .[1:2]
    ndistrib <- rnorm(3,mean = mean(deciles), sd = 1)
}                  