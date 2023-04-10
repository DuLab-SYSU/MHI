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

#################################################################
### Network-based drug repurposing (drug to DEPs after infection) 
### ATC code (2022-12-07)
code_level4 <- read_xlsx('../data/DrugBank/drug_ATC_codes_level4.xlsx')
code_level4 <- rename(code_level4, drug = name)
code_level4$level_4_category <- factor(code_level4$level_4_category,
                                       levels = c(sort(unique(code_level4$level_4_category)[c(2:9,11:15)]), 
                                                  'Various', 'Others'))
code_level4 <- code_level4 %>% 
  arrange(level_4_category)

### for all FDA-approved drugs (2022-12-20)
## dataset from Wang 2022
dir.create("../result/REMAP/drug_proximity_plasma/", recursive = TRUE)

# VTP
vtp_drug <- data.frame()
for (f in list.files(path = "../intermediate/REMAP/drug_proximity_plasma/VTP/")) {
  file_name <- paste0("../intermediate/REMAP/drug_proximity_plasma/VTP/", f)
  drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  vtp_drug <- rbind(vtp_drug, drug_d)
}
write.table(vtp_drug, 
            '../result/REMAP/drug_proximity_plasma/VTP_distance.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

# IN
in_drug <- data.frame()
for (f in list.files(path = "../intermediate/REMAP/drug_proximity_plasma/IN/all/")) {
  file_name <- paste0("../intermediate/REMAP/drug_proximity_plasma/IN/all/", f)
  drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  in_drug <- rbind(in_drug, drug_d)
}
write.table(in_drug, 
            paste0('../result/REMAP/drug_proximity_plasma/IN_distance.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE)

# DEP
dep_drug <- data.frame()
for (f in list.files(path = "../intermediate/REMAP/drug_proximity_plasma/DEP/")) {
  file_name <- paste0("../intermediate/REMAP/drug_proximity_plasma/DEP/", f)
  drug_d <- read.table(file_name, sep = '\t', header = TRUE, quote = '')
  dep_drug <- rbind(dep_drug, drug_d)
}
write.table(dep_drug, 
            paste0('../result/REMAP/drug_proximity_plasma/DEP_distance.txt'),
            quote = FALSE, sep = "\t", row.names = FALSE)

### Candidate drugs (Z score < -1.6)
drug_candidate <- data.frame()
for (f in list.files(path="../result/REMAP/drug_proximity_plasma/")) {
  ca <- str_split(f, "_")[[1]][1]
  
  drugs <- read.table(paste0('../result/REMAP/drug_proximity_plasma/', f),
                      header = TRUE, sep = '\t', quote = '')
  drugs2 <- drugs %>%
    mutate(z.score = as.numeric(z.score)) %>%
    filter(z.score < -1.6) %>%
    mutate(category = ca)
  
  drug_candidate <- rbind(drug_candidate, drugs2)
}

write_xlsx(list(drug_proximity = drug_candidate),
           '../result/REMAP/drug_proximity_plasma/drug_proximity.xlsx')

drug_candidate2 <- aggregate(drug_candidate[6], drug_candidate[1],
                             FUN = function(X) paste(unique(X), collapse=";"))
write.table(drug_candidate2, 
            '../result/REMAP/drug_proximity_plasma/drug_candidate.txt',
            quote = FALSE, sep = "\t", row.names = FALSE)

## Overlap between drugs to those three gene set (2022-12-21)

drug_candidate <- read_xlsx('../result/REMAP/drug_proximity_plasma/drug_proximity.xlsx')

pdf(file = "../result/Fig_5A.pdf")
tid_d_plot <- Venn(list("VTP" = drug_candidate[drug_candidate$category=='VTP',]$drug,
                        "IN" = drug_candidate[drug_candidate$category=='IN',]$drug,
                        "DEP" = drug_candidate[drug_candidate$category=='DEP',]$drug))
plot(tid_d_plot, doWeight = T)
dev.off()

## Circos plot shows a global view of potential drugs across three gene set (2022-12-21)
drug_candidate <- read_xlsx('../result/REMAP/drug_proximity_plasma/drug_proximity.xlsx')

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

pdf('../result/REMAP/drug_proximity_plasma/drug_circos.pdf',
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
remaps_drugs <- read.table('../result/REMAP/REMAP_pred30256.txt',
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

pdf('../result/REMAP/drug_proximity_plasma/drug_category_signif.pdf',
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

