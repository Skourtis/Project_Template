
#BiocManager::install("DEqMS")
library(DEqMS)
url2 <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD000279/proteomebenchmark.zip"
download.file(url2, destfile = "./PXD000279.zip",method = "auto")
unzip("PXD000279.zip")

df.prot = read.table("proteinGroups.txt",header=T,sep="\t",stringsAsFactors = F,
                     comment.char = "",quote ="")

# remove decoy matches and matches to contaminant
df.prot = df.prot[!df.prot$Reverse=="+",]
df.prot = df.prot[!df.prot$Contaminant=="+",]

# Extract columns of LFQ intensites
df.LFQ = df.prot[,89:94]
df.LFQ[df.LFQ==0] <- NA

rownames(df.LFQ) = df.prot$Majority.protein.IDs
df.LFQ$na_count_H = apply(df.LFQ,1,function(x) sum(is.na(x[1:3])))
df.LFQ$na_count_L = apply(df.LFQ,1,function(x) sum(is.na(x[4:6])))
# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_H<2 & df.LFQ$na_count_L<2,1:6]

df.LFQ.filter_var <- df.LFQ.filter %>% rownames_to_column("Uniprot") %>%
    pivot_longer(-Uniprot, names_to = "Experiment", values_to = "Abundance") %>%
    mutate(Abundance = log2(Abundance),
           Experiment = case_when(
               str_detect(Experiment,"H") ~ "Heavy",
               str_detect(Experiment,"L") ~ "Light")) %>%
    na.omit()%>%
    group_by(Uniprot, Experiment) %>%
    summarise(mean_prot = mean(Abundance, na.rm = T),
              variance = var(Abundance, na.rm = T)) 

var_ctr <- ggplot(df.LFQ.filter_var,aes(x = Experiment, y = log2(variance)))+
    geom_boxplot()+
    geom_jitter(alpha = 0.01)

mean_ctr <- ggplot(df.LFQ.filter_var,aes(x = Experiment, y = mean_prot))+
    geom_boxplot()+
    geom_jitter(alpha = 0.01)
library(matrixStats)
# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count.table = data.frame(count = rowMins(as.matrix(df.prot[,19:24])),
                             row.names = df.prot$Majority.protein.IDs)
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count.table$count = pep.count.table$count+1

protein.matrix = log2(as.matrix(df.LFQ.filter))

class = as.factor(c("H","H","H","L","L","L"))
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(protein.matrix,design = design)
cont <- makeContrasts(classH-classL, levels = design)
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, you should troubleshoot the error first
min(fit3$count)
## [1] 1
fit4_ctr = spectraCounteBayes(fit3)
VarianceBoxplot(fit4_ctr, n=30, main = "Label-free dataset PXD000279",
                xlab="peptide count + 1")

df.prot = read.table(here::here("Datasets","Raw","2020MQ044txt","txt","proteinGroups.txt"),header=T,sep="\t",stringsAsFactors = F,
                     comment.char = "",quote ="")

# remove decoy matches and matches to contaminant
df.prot = df.prot[!df.prot$Reverse=="+",]
df.prot = df.prot[!df.prot$`Potential.contaminant`=="+",]

# Extract columns of LFQ intensites
df.LFQ = df.prot %>% select(contains("LFQ"))
df.LFQ[df.LFQ==0] <- NA

rownames(df.LFQ) = df.prot$Majority.protein.IDs
df.LFQ$na_count_H = apply(df.LFQ,1,function(x) sum(is.na(x[1:3])))
df.LFQ$na_count_L = apply(df.LFQ,1,function(x) sum(is.na(x[4:6])))
df.LFQ$na_count_T = apply(df.LFQ,1,function(x) sum(is.na(x[7:9])))
# Filter protein table. DEqMS require minimum two values for each group.
df.LFQ.filter = df.LFQ[df.LFQ$na_count_H<2 & df.LFQ$na_count_L<2 & df.LFQ$na_count_T<2,c(1:9)]
df.LFQ.filter_var_etna <- df.LFQ.filter %>% rownames_to_column("Uniprot") %>%
    pivot_longer(-Uniprot, names_to = "Experiment", values_to = "Abundance") %>%
    mutate(Abundance = log2(Abundance),
           Experiment = case_when(
               str_detect(Experiment,"1|2|3") ~ "Control",
               str_detect(Experiment,"4|5|6") ~ "sh_1",
               str_detect(Experiment,"7|8|9") ~ "sh_2")) %>%
     na.omit()%>%
    group_by(Uniprot, Experiment) %>%
    summarise(mean_prot = mean(Abundance, na.rm = T),
              variance = var(Abundance, na.rm = T)) 
var_etna <-     ggplot(df.LFQ.filter_var_etna,aes(x = Experiment, y = variance))+
        geom_boxplot()+
        geom_jitter(alpha = 0.01)

mean_etna <-     ggplot(df.LFQ.filter_var_etna,aes(x = Experiment, y = mean_prot))+
        geom_boxplot()+
        geom_jitter(alpha = 0.01)
    
library(matrixStats)
# we use minimum peptide count among six samples
# count unique+razor peptides used for quantification
pep.count.table = data.frame(count = rowMins(as.matrix(df.prot[,c(20:22,26:28)])),
                             row.names = df.prot$Majority.protein.IDs)
# Minimum peptide count of some proteins can be 0
# add pseudocount 1 to all proteins
pep.count.table$count = pep.count.table$count+1

protein.matrix = log2(as.matrix(df.LFQ.filter[,1:6]))

class = as.factor(c("H","H","H","L","L","L"))
design = model.matrix(~0+class) # fitting without intercept

fit1 = lmFit(protein.matrix,design = design)
cont <- makeContrasts(classH-classL, levels = design)
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, you should troubleshoot the error first
min(fit3$count)
fit4_etn = spectraCounteBayes(fit3)
VarianceBoxplot(fit4_etn, n=30, main = "Ctrl vs RNF",
                xlab="peptide count + 1")

DEqMS.results = outputResult(fit4,coef_col = 1)
# Add Gene names to the data frame
rownames(df.prot) = df.prot$Majority.protein.IDs
DEqMS.results$Gene.name = df.prot[DEqMS.results$gene,]$Gene.names
head(DEqMS.results %>% arrange(adj.P.Val))


