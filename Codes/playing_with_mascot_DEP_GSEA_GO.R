Mascot <- readxl::read_excel(here::here('Datasets','Raw','2020MQ044-pdresults-table.xlsx'), col_names = F) %>%
    purrr::set_names(.[2,] %>% 
                         unlist()) %>% 
    .[-c(1:2),]
Mascot <- Mascot %>% .[,str_detect(names(.),"Accession|_")] %>%
    .[,-c(2:10)] %>% column_to_rownames('Accession')
Mascot <- Mascot %>%
    mutate(across(where(is.character),as.numeric)) %>%
    na.omit()

#Mascot[is.na(Mascot)] <- 0
install.packages('pheatmap')

# load package
library(pheatmap)
pheatmap(Mascot %>% subset(rownames(.) %in% DEP))

Mascot <- readxl::read_excel(here::here('Datasets','Raw','2020MQ044-pdresults-table.xlsx'), col_names = F) %>%
    purrr::set_names(.[2,] %>% 
                         unlist()) %>% 
    .[-c(1:2),]
Mascot <- Mascot %>% .[,str_detect(names(.),"Accession|adjusted")] %>%
     column_to_rownames('Accession') %>% set_names(c('shR_c','shP_c'))
Mascot <- Mascot %>%
    mutate(across(where(is.character),as.numeric))
DEP <- Mascot %>% 
    subset(shR_c< 0.05 | shP_c< 0.05 ) %>% rownames()
clusterProfiler::gseGO(geneList     = Mascot %>% pull('shR_c') %>% set_names(rownames(Mascot)) %>% sort() %>% .[1:100],
                       OrgDb        = org.Mm.eg.db,
                       ont          = "BP",
                       nPerm        = 1000,
                       minGSSize    = 100,
                       keyType = "UNIPROT",
                       maxGSSize    = 700,
                       pvalueCutoff = 0.05,
                       verbose      = FALSE)
ego <- clusterProfiler::enrichGO(
    gene  = Mascot %>% arrange(shP_c) %>% .[1:100,] %>% rownames() ,
                universe      = rownames(Mascot),
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                keyType = "UNIPROT",
                readable      = TRUE)

barplot(ego, showCategory=20)
clusterProfiler::cnetplot(ego, foldChange= Mascot %>% pull('shR_c') %>% set_names(rownames(Mascot)))
