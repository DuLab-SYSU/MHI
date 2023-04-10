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

### Network-based drug and side effect proximity (REMAP) (2022-12-08)
candi_drug_pred <- read.table('../result/REMAP/drug_proximity/Monkeypox_virus_distance.txt',
                              header = TRUE,
                              sep = '\t',
                              quote = '')
candi_drug_pred <- candi_drug_pred %>%
  inner_join(target_num_pred, by = 'drug') %>%
  filter(z.score < -1.6) %>% 
  arrange(distance, z.score)
candi_drug_pred <- candi_drug_pred %>%  
  mutate(rank = 1:nrow(candi_drug_pred),
         data_set = 'REMAP')

se_proximity_pred <- data.frame()
for (f in list.files(path = "../intermediate/REMAP/drug_sideEffect_proximity/")){
  drug <- read.table(paste0("../intermediate/REMAP/drug_sideEffect_proximity/", f),
                          header = TRUE, sep = '\t', quote = '')
  drug_dist <- drug_dist %>% mutate(Virus = v)
  drug <- drug %>% filter(z.score < -1.6)
  
  se_proximity_pred <- rbind(se_proximity_pred, drug)
}

se_pred_num <- se_proximity_pred %>%
  filter(z.score < -1.6) %>%
  select(drug, side_effect) %>%
  group_by(drug) %>%
  summarise(se_num = n())

## comparison of side-effects of candidate drugs between DrugBank- and REMAP-based datasets (2022-12-08)
se_pred_num <- candi_drug_pred %>%
  select(drug, data_set, overlap) %>%
  left_join(se_pred_num, by = 'drug')

se_proximity <- read.table('../intermediate/DrugBank/drug_sideEffect_proximity/MPXV_277candidate_drugs_sideEffect_proximity.txt',
                           header = TRUE,
                           sep = '\t',
                           quote = '')
se_num <- se_proximity %>%
  filter(z.score < -1.6) %>%
  select(drug, side_effect) %>%
  group_by(drug) %>%
  summarise(se_num = n())
se_num <- candi_drug %>% 
  filter(data_set == 'DrugBank') %>% 
  select(drug, data_set, overlap) %>%
  left_join(se_num, by = 'drug') %>%
  mutate_all(~replace(., is.na(.), 0))

se_num <- rbind(se_num, se_pred_num)
se_num$data_set <- factor(se_num$data_set, levels = c('REMAP', 'DrugBank'))

se_overlap <- se_num %>%
  filter(overlap == 'yes')

library(ggridges)
library(viridis)
#Colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
p_se_distr <- ggplot(se_num, aes(se_num, data_set, fill = data_set)) +
  geom_density_ridges_gradient(alpha = .85) +
  scale_x_continuous(expand = c(0.01, 0), breaks = seq(0, 1200, 400)) +
  scale_y_discrete(expand = c(0.1, 0)) +
  scale_fill_manual(values = c('#E74C3C', '#3498DA')) +
  #scale_fill_viridis(name = "Domain_num", option = "C") +
  #scale_fill_gradientn(colours=Colormap) +
  labs(x = "Side-effect Number", y = '', fill = '') +
  theme_ridges(font_size = 13, grid = FALSE) +
  mytheme +
  theme(axis.title.y = element_blank(),
        legend.position = "none")

ggsave(p_se_distr,
       filename = "../result/Fig_4E.pdf",
       width = 9,
       height = 3,
       units = c("cm"))


library(ggpubr)
p_se_overlap <- ggplot(se_overlap, aes(data_set, se_num)) +
  geom_boxplot(aes(fill = data_set), alpha = .7) +
  scale_fill_manual(values = c('#E74C3C', '#3498DA')) +
  geom_line(aes(group = drug), color = '#2D3D50', lwd = 0.5, alpha = .5) +
  geom_point(size = 3, color = '#2D3D50', shape = 21) +
  labs(x = '', y = 'Side-effect number', fill = '') +
  theme(legend.position = 'top') +
  mytheme +
  theme(legend.position = 'top') +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, 
                     comparisons = list(unique(as.character(se_overlap$data_set))))

# Wilcoxon paired-samples signed rank test
se_num2 <- se_num %>%
  filter(overlap == 'yes') %>%
  select(drug, data_set, se_num) %>%
  pivot_wider(names_from = data_set, values_from = se_num)

wilcox.test(se_num2$REMAP, se_num2$DrugBank, paired = TRUE, alternative = 'greater')


ggsave(p_se_overlap,
       filename = "../result/Fig_4F.pdf",
       width = 3.6,
       height = 5)


### The rank of the candidate drugs (2023-01-11)
drug_proximity <- data.frame()
drug_proximity_signif <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  drug_dist <- read.table(paste0("../result/REMAP/drug_proximity/", v, "_distance.txt"),
                          header = TRUE, sep = '\t', quote = '')
  drug_dist <- drug_dist %>% mutate(Virus = v)
  drug_dist2 <- drug_dist %>% 
    mutate(z.score = as.numeric(z.score)) %>%
    filter(z.score < -1.6)
  
  drug_proximity <- rbind(drug_proximity, drug_dist)
  drug_proximity_signif <- rbind(drug_proximity_signif, drug_dist2)
}

drug_virus_num <- drug_proximity_signif %>%
  group_by(drug) %>%
  summarise(Virus_num = n())

candi_drug_rank <- candi_drug_pred %>%
  select(drug, distance, z.score) %>%
  inner_join(se_pred_num, by = 'drug') %>%
  inner_join(drug_virus_num, by = 'drug') %>% 
  arrange(se_num, Virus_num, z.score, distance, desc(overlap)) %>%
  mutate(Rank = 1:268)

se_num_cutoff <- quantile(candi_drug_rank$se_num, 0.2)
z_cutoff <- quantile(candi_drug_rank$z.score, 0.2)
candi_drug_rank2 <- candi_drug_rank %>%
  mutate(label = ifelse(z.score <= z_cutoff & se_num <= se_num_cutoff, drug, ''),
         color = ifelse(z.score <= z_cutoff & se_num <= se_num_cutoff, '#3498DA', '#BDC3C8'))

write_xlsx(list(drug_proximity = drug_proximity, 
                drug_proximity_zscore = drug_proximity_z,
                candidate_drug_rank = candi_drug_rank2),
           '../result/REMAP/drug_proximity.xlsx')

library(ggrepel)
# Figure 4G
ggplot(candi_drug_rank2, aes(z.score, se_num)) +
  geom_point(aes(size = Virus_num, color = color)) +
  scale_color_manual(values = c('#3498DA', '#BDC3C8')) +
  xlim(-1.631009, -27) +
  ylim(1060, 320) +
  scale_x_break(c(-26, -8), scales = 5) +
  labs(x = 'Z score', y = 'Side-effect Number') +
  geom_text_repel(candi_drug_rank2, mapping = aes(z.score, se_num, label = label), 
                  size = 3, max.overlaps = 200) +
  geom_vline(xintercept = z_cutoff, 
             linetype = 5, color = "#E74C3C") +
  geom_hline(yintercept = se_num_cutoff, 
             linetype = 5, color = "#E74C3C") +
  mytheme


#################################################################
### Network-based drug repurposing (drug to DEPs after infection) 
### for all FDA-approved drugs (2022-12-20)
## dataset from Wang 2022
dir.create("../result/REMAP0.8/drug_proximity_plasma/", recursive = TRUE)

# VTP
vtp_drug <- data.frame()
for (f in list.files(path = "../intermediate/REMAP0.8/drug_proximity_plasma/VTG/")) {
  file_name <- paste0("../intermediate/REMAP0.8/drug_proximity_plasma/VTG/", f)
  drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  vtp_drug <- rbind(vtp_drug, drug_d)
}
write.table(vtp_drug, 
            '../result/REMAP0.8/drug_proximity_plasma/VTP_distance.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

# IN
in_drug <- data.frame()
for (f in list.files(path = "../intermediate/REMAP0.8/drug_proximity_plasma/IN/all/")) {
  file_name <- paste0("../intermediate/REMAP0.8/drug_proximity_plasma/IN/all/", f)
  drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  in_drug <- rbind(in_drug, drug_d)
}
write.table(in_drug, 
            paste0('../result/REMAP0.8/drug_proximity_plasma/IN_distance.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE)

# DEP
dep_drug <- data.frame()
for (f in list.files(path = "../intermediate/REMAP0.8/drug_proximity_plasma/DEP/")) {
  file_name <- paste0("../intermediate/REMAP0.8/drug_proximity_plasma/DEP/", f)
  drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  dep_drug <- rbind(dep_drug, drug_d)
}
write.table(dep_drug, 
            paste0('../result/REMAP0.8/drug_proximity_plasma/DEP_distance.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE)

### Candidate drugs (Z score < -1.6)
drug_candidate <- data.frame()
for (f in list.files(path="../result/REMAP0.8/drug_proximity_plasma/")) {
  ca <- str_split(f, "_")[[1]][1]
  
  drugs <- read.table(paste0('../result/REMAP0.8/drug_proximity_plasma/', f),
                     header = TRUE, sep = '\t', quote = '')
  drugs2 <- drugs %>%
    mutate(z.score = as.numeric(z.score)) %>%
    filter(z.score < -1.6) %>%
    mutate(category = ca)
  
  drug_candidate <- rbind(drug_candidate, drugs2)
}

write_xlsx(list(drug_proximity = drug_candidate),
           '../result/REMAP0.8/drug_proximity_plasma/drug_proximity.xlsx')

drug_candidate2 <- aggregate(drug_candidate[6], drug_candidate[1],
                             FUN = function(X) paste(unique(X), collapse=";"))
write.table(drug_candidate2, 
            '../result/REMAP0.8/drug_proximity_plasma/drug_candidate.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

## Overlap between drugs to those three gene set (2022-12-21)

drug_candidate <- read_xlsx('../result/REMAP0.8/drug_proximity_plasma/drug_proximity.xlsx')

pdf(file = "../result/REMAP0.8/drug_proximity_plasma/drug_proximity.pdf")
tid_d_plot <- Venn(list("VTP" = drug_candidate[drug_candidate$category=='VTP',]$drug,
                        "IN" = drug_candidate[drug_candidate$category=='IN',]$drug,
                        "DEP" = drug_candidate[drug_candidate$category=='DEP',]$drug))
plot(tid_d_plot, doWeight = T)
dev.off()

## Circos plot shows a global view of potential drugs across three gene set (2022-12-21)
drug_candidate <- read_xlsx('../result/REMAP0.8/drug_proximity_plasma/drug_proximity.xlsx')

category_col <- c('#18BC9B', '#2FCC71', '#3498DA', '#9B59B6', '#F1C50F', 
                  '#E67D22', '#E74C3C', '#10A185', '#27AE61', '#2880B9',
                  '#8E45AD', '#F39C12', '#D35301', '#33495E', '#95A6A6')
category_col <- data.frame(level_4_category = levels(code_level4$level_4_category),
                           color = category_col)

code_level4_2 <- code_level4 %>%
  filter(drug %in% unique(drug_candidate$drug)) %>%
  inner_join(category_col, by = 'level_4_category')

drug_3set <- drug_candidate %>%
  select(drug, category) %>%
  mutate(number = 1)

drug_3set2 <- drug_3set %>%
  pivot_wider(names_from = drug, values_from = number) %>%
  column_to_rownames('category') %>%
  mutate_all(~replace(., is.na(.), 0)) 

drug_3set3 <- drug_3set2[rev(rownames(drug_3set2)),code_level4_2$drug]
drug_3set3 <- as.matrix(drug_3set3)


grid.col = NULL
grid.col[colnames(drug_3set3)] <- code_level4_2$color
grid.col[rownames(drug_3set3)] = c('#E74C3C', '#18BC9B', '#F1C50F')

pdf('../result/REMAP0.8/drug_proximity_plasma/drug_circos.pdf',
    width = 5.82,
    height = 8.11)
circos.par(start.degree = 105)
chordDiagram(drug_3set3,
             directional = TRUE,
             annotationTrack = c("name", "grid"),  # 去刻度
             annotationTrackHeight = mm_h(c(3, 3)), # 第二个3调整grid的宽度，
             preAllocateTracks = list(track.height = mm_h(8), track.margin = c(mm_h(4), 0)),
             grid.col = grid.col,
             diffHeight = 0.03,
             transparency = 0.6,
             big.gap = 3,
             small.gap = 0)
dev.off()
circos.clear()


# significant test (2022-12-21)
remaps_drugs <- read.table('../result/REMAP0.8/REMAP_pred30256.txt',
                             header = TRUE,
                             sep = '\t',
                             quote = '')
remaps_drugs <- select(remaps_drugs, drug) %>%
  distinct()

drug_3set3_2 <- drug_3set %>%
  inner_join(code_level4, by = 'drug') %>%
  group_by(category, level_4_category) %>%
  summarise(number = n()) %>%
  pivot_wider(names_from = category, values_from = number)

drug_category_test <- code_level4 %>%
  filter(drug %in% remaps_drugs$drug) %>%
  group_by(level_4_category) %>%
  summarise(Number = n()) %>% 
  inner_join(drug_3set3_2, by = 'level_4_category')

drug_category_test <- drug_category_test %>%
  mutate(DEP_all = sum(drug_category_test$DEP, na.rm = TRUE),
         IN_all = sum(drug_category_test$IN, na.rm = TRUE),
         VTP_all = sum(drug_category_test$VTP, na.rm = TRUE)) %>%
  mutate_all(~replace(., is.na(.), 0))

drug_category_test <- drug_category_test %>%
  mutate(DEP_pvalue = phyper(DEP-1, DEP_all, sum(drug_category_test$Number), Number, lower.tail = F),
         IN_pvalue = phyper(IN-1, IN_all, sum(drug_category_test$Number), Number, lower.tail = F),
         VTP_pvalue = phyper(VTP-1, VTP_all, sum(drug_category_test$Number), Number, lower.tail = F)) %>%
  column_to_rownames('level_4_category')

pdf('../result/REMAP0.8/drug_proximity_plasma/drug_category_signif.pdf',
    width = 4.8,
    height = 6.9)
Heatmap(as.matrix(drug_category_test[,rev(8:10)]),
        name = "P value",
        row_labels = rownames(drug_category_test),
        cluster_columns = FALSE,
        col = colorRamp2(c(min(as.matrix(drug_category_test[,rev(8:10)])), 1), 
                         c("#3398DA", "white")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(as.matrix(drug_category_test[,rev(8:10)])[i, j] < 0.001, "***", 
                           ifelse(as.matrix(drug_category_test[,rev(8:10)])[i, j] < 0.01, "**",
                                  ifelse(as.matrix(drug_category_test[,rev(8:10)])[i, j] < 0.05, "*", ''))),
                    x, y, gp = gpar(fontsize = 10))
        },
        row_names_gp = gpar(col = c('#18BC9B', '#2FCC71', '#3498DA', '#9B59B6', '#F1C50F',
                                    '#E67D22', '#E74C3C', '#10A185', '#27AE61', '#2880B9',
                                    '#8E45AD', '#F39C12', '#D35301', '#33495E', '#95A6A6')),
        column_labels = c('VTP', 'IN', 'DEP'))

dev.off()

pdf('../result/REMAP0.8/drug_proximity_plasma/drug_category_signif2.pdf',
    width = 6.9,
    height = 4.2)

Heatmap(t(as.matrix(drug_category_test[,rev(8:10)])),
        name = "P value",
        cluster_rows = FALSE,
        col = colorRamp2(c(min(as.matrix(drug_category_test[,rev(8:10)])), 1), 
                         c("#3398DA", "white")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(t(as.matrix(drug_category_test[,rev(8:10)]))[i,j] < 0.001, "***",
                           ifelse(t(as.matrix(drug_category_test[,rev(8:10)]))[i,j] < 0.01, "**",
                                  ifelse(t(as.matrix(drug_category_test[,rev(8:10)]))[i,j] < 0.05, "*", ''))),
                    x, y, gp = gpar(fontsize = 10))
        },
        column_names_gp = gpar(col = c('#18BC9B', '#2FCC71', '#3498DA', '#9B59B6', '#F1C50F',
                                       '#E67D22', '#E74C3C', '#10A185', '#27AE61', '#2880B9',
                                       '#8E45AD', '#F39C12', '#D35301', '#33495E', '#95A6A6')),
        row_labels = c('VTP', 'IN', 'DEP'))

dev.off()

### Final ranked drug list
## Number of side-effects for the candidate drugs (2022-12-26)
drug_plasma_se <- data.frame()
for (f in list.files(path = "../intermediate/REMAP0.8/drug_sideEffect_proximity_plasma/")) {
  file_name <- paste0("../intermediate/REMAP0.8/drug_sideEffect_proximity_plasma/", f)
  drug_se <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  drug_plasma_se <- rbind(drug_plasma_se, drug_se)
}

drug_plasma_se <- drug_plasma_se %>%
  filter(z.score < -1.6)
write_xlsx(list(drug_plasma_se = drug_plasma_se),
           '../result/REMAP0.8/drug_proximity_plasma/drug_plasma_se.xlsx')
  

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
       filename = "../result/REMAP0.8/drug_proximity_plasma/drug_candidate_3sets.pdf",
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

pdf('../result/REMAP0.8/drug_proximity_plasma/drug_category_filtered_signif.pdf',
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

### Network-based proximity between drug and pathways (2023-02-08)
drug_pathway_proximity <- read.table('../result/REMAP0.8/drug_pathway_proximity_plasma/MPXV_candidate_drugs_pathway_proximity.txt',
                                     header = TRUE, sep = '\t', quote = '')
drug_pathway_proximity <- drug_pathway_proximity %>%
  filter(z.score < -1.6) %>% 
  select(drug, pathway)

drug_pathway_proximity <- drug_candidate %>%
  select(drug, distance, z.score, category) %>%
  inner_join(drug_pathway_proximity, by = 'drug')

p_pathway <- read_xlsx('../result/MPXV_VTG&IN&Wang_DEP_reactome_pathway_ALL2.xlsx')
p_pathway <- p_pathway %>%
  select(Category, Description, Pathway_category) %>% 
  mutate(Category = str_replace(Category, 'VTG', 'VTP'))
colnames(p_pathway) <- c('category', 'pathway', 'pathway_category')

drug_pathway_proximity2 <- drug_pathway_proximity %>%
  inner_join(p_pathway, by = c('category', 'pathway')) %>%
  select(drug, pathway_category) %>%
  distinct()
drug_pathway_proximity2 <- aggregate(drug_pathway_proximity2[2], drug_pathway_proximity2[1],
                                     FUN = function(X) paste(unique(X), collapse=";"))


final_candidate <- left_join(drug_plasma_se_num, drug_pathway_proximity2, by = 'drug') %>%
  arrange(desc(protein_set_num), se_num, distance)
final_candidate <- candi_3sets_drug_rank2[,c(1,3,4)] %>%
  inner_join(final_candidate, by = 'drug')

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
           '../result/REMAP0.8/drug_proximity_plasma/final_candidates.xlsx')

