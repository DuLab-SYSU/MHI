library(RColorBrewer)
library(tidyverse)
library(readxl)
library(writexl)

library(ComplexHeatmap)
library(circlize)
library(Vennerable)

mytheme <- theme_classic() +
  theme(panel.grid.minor = element_line(colour=brewer.pal(9,"Pastel1")[9],
                                        linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        plot.margin = unit(rep(0.3, 4),'cm'))

setwd('~/Desktop/Poxviridae-human/src/')

### ATC code (2022-12-07)
code_level4 <- read_xlsx('../data/DrugBank/drug_ATC_codes_level4.xlsx')
code_level4 <- rename(code_level4, drug = name)
code_level4$level_4_category <- factor(code_level4$level_4_category,
                                       levels = c(sort(unique(code_level4$level_4_category)[c(2:9,11:15)]), 
                                                  'Various', 'Others'))
code_level4 <- code_level4 %>% 
  arrange(level_4_category)

### Final ranked drug list
drug_candidate <- read_xlsx('../result/REMAP/drug_proximity_plasma/drug_proximity.xlsx')
drug_candidate2 <- aggregate(drug_candidate[6], drug_candidate[1],
                             FUN = function(X) paste(unique(X), collapse=";"))

## Number of side-effects for the candidate drugs (2022-12-26)
drug_plasma_se <- data.frame()
for (f in list.files(path = "../intermediate/REMAP/drug_sideEffect_proximity_plasma/")) {
  file_name <- paste0("../intermediate/REMAP/drug_sideEffect_proximity_plasma/", f)
  drug_se <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  drug_plasma_se <- rbind(drug_plasma_se, drug_se)
}

drug_plasma_se <- drug_plasma_se %>%
  filter(z.score < -1.6)
write_xlsx(list(drug_plasma_se = drug_plasma_se),
           '../result/REMAP/drug_proximity_plasma/drug_plasma_se.xlsx')


drug_plasma_se_num <- drug_plasma_se %>%
  group_by(drug) %>%
  summarise(se_num = n()) %>%
  right_join(drug_candidate2, by = "drug") %>%
  mutate_all(~replace(., is.na(.), 0))

drug_plasma_se_num <- drug_candidate %>% 
  select(drug, distance) %>%
  group_by(drug) %>% 
  slice_min(distance) %>% 
  distinct() %>%
  inner_join(drug_plasma_se_num, by = 'drug') %>%
  arrange(se_num, distance) %>%
  rename(protein_set = category) %>%
  mutate(protein_set_num = length(str_split(protein_set, ';')[[1]]))

drug_plasma_se_num <- drug_plasma_se_num %>%
  inner_join(code_level4, by = "drug")

## The narrow down of the drug list (2023-02-08)
candi_3sets_drug_rank <- drug_plasma_se %>%
  group_by(drug) %>%
  summarise(se_num = n()) %>%
  right_join(drug_candidate, by = "drug") %>%
  mutate_all(~replace(., is.na(.), 0))

se_num_cutoff2 <- quantile(candi_3sets_drug_rank$se_num, 0.2)
z_cutoff2 <- quantile(candi_3sets_drug_rank$z.score, 0.2)
candi_3sets_drug_rank2 <- candi_3sets_drug_rank %>%
  mutate(label = ifelse(z.score <= z_cutoff2 & se_num <= se_num_cutoff2, drug, ''))

candi_3sets_drug_rank2 <- candi_3sets_drug_rank2 %>% 
  select(drug, se_num, z.score, label) %>% 
  group_by(drug, se_num) %>% 
  slice_min(z.score) %>%
  inner_join(drug_candidate2, by = 'drug') %>%
  mutate(color = ifelse(str_length(label) == 0, '#BDC3C8', category))
candi_3sets_drug_rank2$protein_set_num <- apply(candi_3sets_drug_rank2[,5,drop = FALSE], 1, 
                                                function(x) length(str_split(x, ";")[[1]]))
candi_3sets_drug_rank2 <- candi_3sets_drug_rank2 %>%
  mutate(protein_set_num = as.character(protein_set_num))

p_3sets_drug <- ggplot(candi_3sets_drug_rank2, aes(z.score, se_num)) +
  geom_point(aes(size = protein_set_num, color = color)) +
  scale_size_manual(values = c(2,4,6)) +
  scale_color_manual(values = c('#BDC3C8', '#E8BF22', '#EE7C6F', '#26B394', '#DA4D3C')) +
  xlim(max(candi_3sets_drug_rank2$z.score), min(candi_3sets_drug_rank2$z.score)) +
  ylim(max(candi_3sets_drug_rank2$se_num), min(candi_3sets_drug_rank2$se_num)) +
  labs(x = 'Z score', y = 'Side-effect Number') +
  geom_vline(xintercept = z_cutoff2, 
             linetype = 5, color = "#E74C3C") +
  geom_hline(yintercept = se_num_cutoff2, 
             linetype = 5, color = "#E74C3C") +
  mytheme +
  theme(legend.position = c(0.8, 0.5))
ggsave(p_3sets_drug,
       filename = "../result/Fig_5B.pdf",
       width = 8,
       height = 8,
       units = c("cm"))

## ATC enrichment of filtered drugs (2023-02-08)
drug_3set_filter <- candi_3sets_drug_rank %>%
  filter(z.score <= z_cutoff2 & se_num <= se_num_cutoff2) %>%
  select(drug, category) %>%
  mutate(number = 1)

drug_3set_filter2 <- drug_3set_filter %>%
  inner_join(code_level4, by = 'drug') %>%
  group_by(category, level_4_category) %>%
  summarise(number = n()) %>%
  pivot_wider(names_from = category, values_from = number)

drug_category_filter_test <- code_level4 %>%
  filter(drug %in% remaps_drugs$drug) %>%
  group_by(level_4_category) %>%
  summarise(Number = n()) %>% 
  inner_join(drug_3set_filter2, by = 'level_4_category')

drug_category_filter_test <- drug_category_filter_test %>%
  mutate(DEP_all = sum(drug_category_filter_test$DEP, na.rm = TRUE),
         IN_all = sum(drug_category_filter_test$IN, na.rm = TRUE),
         VTP_all = sum(drug_category_filter_test$VTP, na.rm = TRUE)) %>%
  mutate_all(~replace(., is.na(.), 0))

drug_category_filter_test <- drug_category_filter_test %>%
  mutate(DEP_pvalue = phyper(DEP-1, DEP_all, sum(drug_category_filter_test$Number), Number, lower.tail = F),
         IN_pvalue = phyper(IN-1, IN_all, sum(drug_category_filter_test$Number), Number, lower.tail = F),
         VTP_pvalue = phyper(VTP-1, VTP_all, sum(drug_category_filter_test$Number), Number, lower.tail = F)) %>%
  column_to_rownames('level_4_category')

pdf('../result/Fig_5C.pdf',
    width = 4.8,
    height = 4.4)
Heatmap(t(as.matrix(drug_category_filter_test[,rev(8:10)])),
        name = "P value",
        cluster_rows = FALSE,
        col = colorRamp2(c(min(as.matrix(drug_category_filter_test[,rev(8:10)])), 1), 
                         c("#3398DA", "white")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(t(as.matrix(drug_category_filter_test[,rev(8:10)]))[i,j] < 0.001, "***",
                           ifelse(t(as.matrix(drug_category_filter_test[,rev(8:10)]))[i,j] < 0.01, "**",
                                  ifelse(t(as.matrix(drug_category_filter_test[,rev(8:10)]))[i,j] < 0.05, "*", ''))),
                    x, y, gp = gpar(fontsize = 10))
        },
        row_labels = c('VTP', 'IN', 'DEP'))

dev.off()

### Final drug list (2023-02-08)

final_candidate <- candi_3sets_drug_rank2[,c(1,3,4)] %>%
  inner_join(drug_plasma_se_num, by = 'drug')

three_set_candidate <- final_candidate %>%
  filter(protein_set == 'DEP;IN;VTP')
DEPandIN_candidate <- final_candidate %>%
  filter(protein_set == 'DEP;IN')
DEPandVTP_candidate <- final_candidate %>%
  filter(protein_set == 'DEP;VTP')
INandVTP_candidate <- final_candidate %>%
  filter(protein_set == 'IN;VTP')
VTP_candidate <- final_candidate %>%
  filter(protein_set == 'VTP')
IN_candidate <- final_candidate %>%
  filter(protein_set == 'IN')
DEP_candidate <- final_candidate %>%
  filter(protein_set == 'DEP')

write_xlsx(list(`DEP;IN;VTP` = three_set_candidate,
                `DEP;IN` = DEPandIN_candidate,
                `DEP;VTP` = DEPandVTP_candidate,
                `IN;VTP` = INandVTP_candidate,
                `VTP` = VTP_candidate,
                `IN` = IN_candidate,
                `DEP` = DEP_candidate),
           '../result/final_candidates.xlsx')
