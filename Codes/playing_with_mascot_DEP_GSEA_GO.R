pacman::p_load(tidyverse,DEqMS,matrixStats)
Mascot <- readxl::read_excel(here::here('Datasets','Raw','2020MQ044-pdresults-table.xlsx'), col_names = F) %>%
    purrr::set_names(.[2,] %>% 
                         unlist()) %>% 
    .[-c(1:2),]

####testing 0 pvalues
Mascot <- Mascot %>% select(contains("pval"))



adjusted <- c(rep(0,100),runif(1000,0,0.02),runif(5000,0,1)) %>% #,
    p.adjust(method = "BH")

adjusted %>% .[.<0.05] %>% length()
adjusted %>% hist(xlim = c(0,1) )

###
Mascot <- Mascot %>% .[,str_detect(names(.),"Accession|_")] %>%
    .[,-c(2:10)] %>% column_to_rownames('Accession')
row_names <- rownames(Mascot)
Mascot <- Mascot %>% 
    mutate(across(where(is.character),as.numeric)) 
rownames(Mascot)  <- row_names   
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

Imp_mascot <- Proteomics_imputation(Data = Mascot,N_condition = 3, N_rep = 3)
Imp_mascot %>% is.na() %>% table()
Mascot %>% is.na() %>% table()


df.LFQ.filter_var_etna <- Imp_mascot %>% as.data.frame() %>% rownames_to_column("Uniprot") %>%
    pivot_longer(-Uniprot, names_to = "Experiment", values_to = "Abundance") %>%
    mutate(Abundance = log2(Abundance),
           Experiment = case_when(
               str_detect(Experiment,"_01|_02|_03") ~ "Control",
               str_detect(Experiment,"_04|_05|_06") ~ "sh_1",
               str_detect(Experiment,"_07|_08|_09") ~ "sh_2")) %>%
    na.omit()%>%
    group_by(Uniprot, Experiment) %>%
    summarise(mean_prot = mean(Abundance, na.rm = T),
              variance = var(Abundance, na.rm = T)) 

var_etna <-     ggplot(df.LFQ.filter_var_etna,aes(x = Experiment, y = log2(variance)))+
    geom_boxplot()+
    geom_jitter(alpha = 0.01)

mean_etna <-     ggplot(df.LFQ.filter_var_etna,aes(x = Experiment, y = mean_prot))+
    geom_boxplot()+
    geom_jitter(alpha = 0.01)
df.prot = read.table(here::here("Datasets","Raw","2020MQ044txt","txt","proteinGroups.txt"),header=T,sep="\t",stringsAsFactors = F,
                     comment.char = "",quote ="")

# remove decoy matches and matches to contaminant
df.prot = df.prot[!df.prot$Reverse=="+",]
df.prot = df.prot[!df.prot$`Potential.contaminant`=="+",]

# Extract columns of LFQ intensites
df.LFQ = df.prot %>% select(contains("LFQ"))
df.LFQ[df.LFQ==0] <- NA
rownames(df.LFQ) = df.prot$Majority.protein.IDs
df.LFQ <- df.LFQ %>% Proteomics_imputation(N_condition = 3, N_rep = 3) %>%
    as.data.frame()
#df.LFQ_imp
pep.count.table = data.frame(count = rowMins(as.matrix(df.prot[,c(23:25,26:29)])),
                             row.names = df.prot$Majority.protein.IDs)
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count.table$count = pep.count.table$count+1

#df.LFQ <- Imp_mascot %>% as.data.frame()
df.LFQ$na_count_C = apply(df.LFQ,1,function(x) sum(is.na(x[1:3])))
df.LFQ$na_count_R = apply(df.LFQ,1,function(x) sum(is.na(x[4:6])))
df.LFQ$na_count_Tp = apply(df.LFQ,1,function(x) sum(is.na(x[7:9])))

df.LFQ.filter = df.LFQ[df.LFQ$na_count_C<2 & df.LFQ$na_count_R<2 & df.LFQ$na_count_Tp<2,c(1:9)] 
protein.matrix_CvsR = log2(as.matrix(df.LFQ.filter[,c(1:3,7:9)]))
class = as.factor(c('H','H','H','L','L','L'))
design =model.matrix(~0+class)

fit1CvsR = lmFit(protein.matrix_CvsR,design = design)
cont <- makeContrasts(classH-classL, levels = design)
fit2CvsR = contrasts.fit(fit1CvsR,contrasts = cont)
fit3CvsR <- eBayes(fit2CvsR)

fit3CvsR$count = pep.count.table[rownames(fit3CvsR$coefficients),"count"]

min(fit3CvsR$count)
fit4 = spectraCounteBayes(fit3CvsR)

VarianceBoxplot(fit4, n=20, main = "Label-free dataset PXD000279",
                xlab="peptide count + 1")

DEqMS.results = outputResult(fit4,coef_col = 1)
# Add Gene names to the data frame
rownames(df.prot) = df.prot$Majority.protein.IDs
DEqMS.results$Gene.name = df.prot[DEqMS.results$gene,]$Gene.names
head(DEqMS.results)
