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
viruses <- c("Cowpox_virus", "Ectromelia_virus-ERPV", "Horsepox_virus", 
             "Molluscum_contagiosum_virus-subtype-1", "Monkeypox_virus", "Monkeypox_virus-Zaire-96-I-16",
             "Orf_virus", "Vaccinia_virus", "Vaccinia_virus-Ankara", 
             "Variola_virus", "Yaba-like_disease_virus", "Yaba_monkey_tumor_virus")
abbre <- c("CPXV", "ECTV ERPV", "HPXV", "MCV1", "MPXV", "MPXV Zaire-96-I-16",
           "ORFV", "VACV", "VACV Ankara", "VARV", "YLDV", "YMTV")

virus_abbre <- data.frame(Virus = viruses,
                          Abbreviation = abbre)

library(igraph)
ppi <- read.table('../result/integrated_ppi_largest_component.txt',
                  header = FALSE,
                  sep = '\t',
                  quote = '')
ppin <- graph_from_data_frame(ppi, directed = FALSE)
nodes <- get.vertex.attribute(ppin)
nodes <- nodes$name

#########################################################
### Network-based drug repurposing (drug to viral module)
### ATC code (2022-12-07)
code_level4 <- read_xlsx('../data/DrugBank/drug_ATC_codes_level4.xlsx')
code_level4 <- rename(code_level4, drug = name)
code_level4$level_4_category <- factor(code_level4$level_4_category,
                                       levels = c(sort(unique(code_level4$level_4_category)[c(2:9,11:15)]), 
                                                  'Various', 'Others'))
code_level4 <- code_level4 %>% 
  arrange(level_4_category)

### for all FDA-approved drugs (2022-12-13)
dir.create("../result/REMAP0.8/drug_proximity/", recursive = TRUE)
for (dire in list.dirs(path = "../intermediate/REMAP/drug_proximity/", recursive = FALSE)) {
  virus <- stringr::str_split(dire, "/")[[1]][6]
  
  drug_dist <- data.frame()
  for (d in list.files(path = dire)){
    file_name <- paste0(dire, "/", d)
    drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
    drug_dist <- rbind(drug_dist, drug_d)
  }
  write.table(drug_dist, 
              paste0('../result/REMAP/drug_proximity/', virus, "_distance.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}

### comparison target number of candidate drugs between DrugBank- and REMAP-based datasets (2022-12-09) 
dti <- read.table('../data/REMAP/drug_targets_approved.tsv',
                  header = TRUE,
                  sep = '\t',
                  quote = '') %>%
  rename(protein = gene) %>%
  filter(protein %in% nodes)
target_num <- dti %>%
  group_by(drug) %>%
  summarise(target_num = n())

dti_pred <- read.table('../result/REMAP/REMAP_pred30256.txt',
                       header = TRUE,
                       sep = '\t',
                       quote = '')
target_num_pred <- dti_pred %>%
  filter(protein %in% nodes) %>%
  group_by(drug) %>%
  summarise(target_num = n())

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


candi_drug <- read.table('../result/DrugBank/drug_proximity/Monkeypox_virus_distance.txt',
                         header = TRUE,
                         sep = '\t',
                         quote = '')
candi_drug <- candi_drug %>%
  inner_join(target_num, by = 'drug') %>%
  filter(z.score < -1.6) %>% 
  arrange(distance, z.score)
candi_drug <- candi_drug %>%  
  mutate(rank = 1:nrow(candi_drug),
         data_set = 'DrugBank') %>%
  mutate(overlap = ifelse(drug %in% candi_drug_pred$drug, "yes", "no"))

candi_drug_pred <- candi_drug_pred %>%
  mutate(overlap = ifelse(drug %in% candi_drug$drug, "yes", "no"))

candi_drug <- rbind(candi_drug, candi_drug_pred)
candi_drug$overlap <- factor(candi_drug$overlap, levels = c('yes', 'no'))
candi_drug$data_set <- factor(candi_drug$data_set, levels = c('REMAP', 'DrugBank'))

candi_drug_overlap <- candi_drug %>%
  filter(overlap == 'yes')

# Wilcoxon paired-samples signed rank test
target_num2 <- candi_drug_overlap %>%
  select(drug, data_set, target_num) %>%
  pivot_wider(names_from = data_set, values_from = target_num)

wilcox.test(target_num2$REMAP, target_num2$DrugBank, paired = TRUE, alternative = 'greater')

pdf(file = "../result/Fig_4D.pdf")
drug_comp_plot <- Venn(list("REMAP" = candi_drug[candi_drug$data_set == 'REMAP',1],
                            "DrugBank" = candi_drug[candi_drug$data_set == 'DrugBank',1]))
plot(drug_comp_plot, doWeight = T)
dev.off()

fisher.test(matrix(c(71, 197, 206, 2298-474), 
                   nrow = 2), alternative = "two.sided")

# candidate drug list (2022-12-15)
candi_drug2 <- aggregate(candi_drug[8], candi_drug[1],
                         FUN = function(X) paste(unique(X), collapse=";"))
write.table(candi_drug2, 
            '../result/REMAP/candidate_drugs.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)


### Number of drugs prioritized across poxviruses and number of viruses for each drugs (2022-12-13)
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

drug_proximity <- inner_join(drug_proximity, virus_abbre, by = 'Virus') %>%
  inner_join(code_level4, by = 'drug')
drug_proximity_signif <- inner_join(drug_proximity_signif, virus_abbre, by = 'Virus') %>%
  inner_join(code_level4, by = 'drug')

drug_num <- drug_proximity_signif %>%
  group_by(Abbreviation) %>%
  summarise(Number = n()) %>%
  arrange(desc(Number))

virus_num <- drug_proximity_signif %>%
  group_by(drug) %>%
  summarise(Virus_num = n()) %>%
  group_by(Virus_num) %>%
  summarise(Drug_num = n())


p_num <- ggplot(drug_num, aes(Abbreviation, Number)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(x = Abbreviation, y = Number - 20, label = Number),
            size = 3, color = "white", angle = 90) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x="",y="Number of drugs") +
  scale_x_discrete(limits = drug_num$Abbreviation) +
  mytheme +
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

p_num2 <- ggplot(virus_num, aes(factor(Virus_num, levels = 1:12), Drug_num)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(x = Virus_num, y = Drug_num - 5, label = Drug_num),
            size = 3, color = "white", angle = 90) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x="Number of poxviruses",y="Number of drugs") +
  mytheme

p_virus_drug_num <- ggpubr::ggarrange(p_num, p_num2,
                                      nrow=1, labels = c("A", "B"),
                                      font.label = list(size = 11))
ggsave(p_virus_drug_num,
       filename = "../result/Fig_4AB.pdf",
       width = 15,
       height = 7,
       units = c("cm"))

### heatmap highlighting network proximity value for 12 poxviruses (2022-12-14)
colors2 <- c('#18BC9B', '#2FCC71', '#3498DA', '#9B59B6', '#95A6A6',
             '#F1C50F', '#E67D22', '#E74C3C', '#33495E', '#BDC3C8')

drug_proximity_wide <- drug_proximity %>%
  dplyr::select(drug, level_4, level_4_category, distance, Abbreviation) %>%
  pivot_wider(names_from = Abbreviation, values_from = distance)

drug_proximity_z <- drug_proximity %>%
  dplyr::select(drug, level_4, level_4_category, z.score, Abbreviation) %>%
  pivot_wider(names_from = Abbreviation, values_from = z.score) %>%
  filter(MPXV < -1.6) %>%
  arrange(MPXV)

drug_proximity_z_long <- drug_proximity %>%
  filter(drug %in% drug_proximity_z$drug) %>%
  mutate(p.value = ifelse(z.score < -1.6, '<0.05', '>=0.05'))

### Pearson correlation coefficients of the proximities for each viral pair (2022-12-14)
pear_cor <- function(proximity_pearson, c1, c2){
  values <- cor.test(as.numeric(as.matrix(drug_proximity_wide[proximity_pearson[c1]])), 
                     as.numeric(as.matrix(drug_proximity_wide[proximity_pearson[c2]])), 
                     method="pearson")
  values$estimate
}

proximity_pearson <- as.data.frame(t(combn(abbre, 2)))
pear_p <- function(proximity_pearson, c1, c2){
  values <- cor.test(as.numeric(as.matrix(drug_proximity_wide[proximity_pearson[c1]])), 
                     as.numeric(as.matrix(drug_proximity_wide[proximity_pearson[c2]])), 
                     method="pearson")
  values$p.value
}

colnames(proximity_pearson) <- c("virus_A", "virus_B")
proximity_pearson$correlation <- apply(proximity_pearson[,1:2], 1, 
                                       pear_cor, c1 = 'virus_A', c2 = 'virus_B')
proximity_pearson$p_value <- apply(proximity_pearson[,1:2], 1, 
                                   pear_p, c1 = 'virus_A', c2 = 'virus_B')

proximity_pearson2 <- proximity_pearson %>%
  dplyr::select(virus_B, virus_A, correlation, p_value) %>%
  dplyr::rename(virus_A = virus_B, virus_B = virus_A) %>%
  rbind(proximity_pearson) %>%
  dplyr::select(-p_value) %>%
  pivot_wider(names_from = virus_B, values_from = correlation) %>%
  arrange(virus_A)
proximity_pearson2[is.na(proximity_pearson2)] <- 1

library(ComplexHeatmap)
library(circlize)
pdf('../result/Fig_4C.pdf',
    width = 10,
    height = 7)
Heatmap(as.matrix(proximity_pearson2[,2:13]),
        name = 'Pearson correlation coefficient',
        row_labels = proximity_pearson2$virus_A,
        col = colorRamp2(c(min(as.matrix(proximity_pearson2[,2:13])), 
                           max(as.matrix(proximity_pearson2[,2:13]))), 
                         c("white", "#E74C3C")))
dev.off()