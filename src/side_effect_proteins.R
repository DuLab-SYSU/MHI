library(tidyverse)

setwd('~/Documents/Computational_Biology/Poxviridae-human/src/')
### drug targets
drug_target <- read.table('../data/DrugBank/drug_targets_approved.txt',
                          header = TRUE,
                          sep = '\t',
                          quote = '') 
drug_target <- select(drug_target, Name, Gene_target, ATC_code) %>%
  rename(name = Name, protein = Gene_target) %>%
  distinct()

### drug and side effect
drug_atc <- read.table('../data/SIDER/drug_atc.tsv',
                        header = FALSE,
                        sep = '\t',
                        quote = '')
colnames(drug_atc) <- c('id', 'ATC_code')

meddra_se <- read.table('../data/SIDER/meddra_all_se.tsv',
                        header = FALSE,
                        sep = '\t',
                        quote = '')
drug_se <- filter(meddra_se, V4 == 'PT') %>%
  select(V1, V6) %>%
  rename(id = V1, side_effect = V6) %>%
  inner_join(drug_atc, by = 'id') %>%
  select(ATC_code, side_effect) %>%
  inner_join(drug_target, by = 'ATC_code') %>%
  select(name, protein, side_effect) %>%
  distinct()

### side-effect proteins
protein_drug_num <- select(drug_se, name, protein) %>%
  distinct() %>%
  group_by(protein) %>%
  summarise(protein_drug_num = n())

se_drug_num <- select(drug_se, name, side_effect) %>%
  distinct() %>%
  group_by(side_effect) %>%
  summarise(se_drug_num = n())

protein_se <- drug_se %>%
  group_by(protein, side_effect) %>%
  summarise(both_num = n()) %>%
  inner_join(protein_drug_num, by = 'protein') %>%
  inner_join(se_drug_num, by = 'side_effect') %>%
  mutate(neither_num = 598 - protein_drug_num - se_drug_num + both_num) %>%
  filter(protein_drug_num >= 5 & se_drug_num >= 5)

fisherTest <- function(x1, x2, x3, x4){
  tmp<- fisher.test(matrix(c(x1, x2-x1, x3-x1, x4), 
                           nrow = 2), alternative = "greater")
  return(tmp$p.value)
}

protein_se <- as.data.frame(protein_se)
protein_se$p_value <- apply(protein_se[,3:6], 1, 
                            function(x) fisher.test(matrix(c(x[1], x[3]-x[1], x[2]-x[1], x[4]), 
                                                           nrow = 2), alternative = "greater")$p.value)

protein_se <- protein_se %>%
  mutate(ADJ_pvalue = p.adjust(p_value, method="BH", n = nrow(protein_se))) %>%
  filter(ADJ_pvalue < 0.2) %>%
  select(side_effect, protein, ADJ_pvalue)

write.table(protein_se, 
            '../data/SIDER/side_effect_proteins.txt',
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
