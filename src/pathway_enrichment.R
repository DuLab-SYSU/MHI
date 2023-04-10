###############################################
# Reactome pathway over-representation analysis

library(clusterProfiler)
library(tidyverse)
library(readxl)
library(writexl)
library(RColorBrewer)
mytheme <- theme(panel.grid.major = element_line(colour=brewer.pal(9,"Pastel1")[9],
                                                 linetype = "longdash"),
                 panel.background = element_rect(fill='transparent', color="#000000"),
                 panel.border=element_rect(fill='transparent', color='black'),
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 9),
                 legend.text = element_text(size = 7),
                 legend.title = element_text(size = 8),
                 legend.background = element_blank())

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

jaccard <- function(x,y){
  jaccard_index <- length(intersect(x, y))/length(union(x, y))
  jaccard_index
}

### VTGs, INs and DEPs (from Wang et al. 2022) (2022-10-28)
plasma_ppi <- read.table('../result/plasma_specific_network.txt',
                         header = FALSE,
                         sep = '\t',
                         quote = '')
plasma_net <- graph_from_data_frame(plasma_ppi, directed = FALSE)
pnodes <- get.vertex.attribute(plasma_net)
pnodes <- pnodes$name

# VTGs in plasma network
MPXV_ppi <- read.table('../data/p-hipster/Monkeypox_virus.txt',
                       header = TRUE,
                       sep = '\t',
                       quote = '')
MPXV_ppi <- filter(MPXV_ppi, Human.Protein %in% nodes)
vtgs <- unique(MPXV_ppi$Human.Protein)
p_vtgs <- vtgs[vtgs %in% pnodes]

# DEPs in plasma network
deps <- read.table('../data/proteome/MPXV_plasma_DEP.txt',
                   header = FALSE,
                   sep = '\t',
                   quote = '')

deps <- filter(deps, V1 %in% pnodes)

## Reactome pathway enrichment
pathway <- read.table('../data/Reactome_pathway/reactome_pathways_from_MsigDB.txt',
                      header = TRUE, sep = '\t', quote = '')
pathway <- filter(pathway, Gene %in% pnodes) %>%
  select(Brief_description, Gene) %>%
  rename(Description = Brief_description)

for (reg in c('all', 'up', 'down')) {
  ins <- read.table(paste0('../intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/intermediate_nodes_', reg, '.txt'),
                    header = FALSE,
                    sep = '\t')
  ins <- ins$V1
  
  if(reg=='all'){
    reg_deps <- deps$V1
  } else {
    reg_deps <- deps[deps$V2==reg,1]
  }
  
  p_net_genes <- data.frame(Gene = c(p_vtgs, ins, reg_deps),
                            Category = c(rep('VTG', length(p_vtgs)), rep('IN', length(ins)), rep('DEP', length(reg_deps))))
  
  p_net_pathway <- data.frame()
  for(ca in unique(p_net_genes$Category)){
    gene_list <- unique(p_net_genes[p_net_genes$Category==ca,1])
    print(length(gene_list))
    if (length(intersect(gene_list, unique(pathway$Gene))) == 0) next
    reactome <- enricher(gene_list,
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   universe = pnodes,
                   qvalueCutoff = 0.01,
                   minGSSize = 1,
                   maxGSSize = 50000,
                   TERM2GENE= pathway)
    reactome2 <- reactome@result
    reactome2 <- reactome2[reactome2$pvalue < 0.01 & reactome2$p.adjust < 0.01,-1]
    reactome2 <- reactome2[,1:7]
    reactome2$Category <- rep(ca, nrow(reactome2))
    
    if(nrow(reactome2)>0){
      reduRow <- c()
      for(n in 1:nrow(reactome2)){
        for(m in 1:nrow(reactome2)){
          if(n != m){
            list1 = strsplit(reactome2$geneID[n], split = '/', fixed = T)[[1]]
            list2 = strsplit(reactome2$geneID[m], split = '/', fixed = T)[[1]]
            ji = jaccard(list1, list2)
            if(ji >= 0.9){
              if(reactome2$p.adjust[n] > reactome2$p.adjust[m]){
                reduRow <- c(reduRow, n)
              } else {
                reduRow <- c(reduRow, m)
              }
            }
          }
        }
      }
      
      if(length(reduRow)>0){
        reactome3 <- reactome2[-unique(reduRow),]
        p_net_pathway <- rbind(p_net_pathway, reactome3)
      }else{
        p_net_pathway <- rbind(p_net_pathway, reactome2)
      }
    }
    print(paste("Category", ca, "finished!"))
  }
  
  p_net_pathway <- p_net_pathway %>% 
    mutate(GeneRatio = paste0("(", GeneRatio, ")"), BgRatio = paste0("(", BgRatio, ")")) %>%
    select(Category, Description,	GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID)
  
  write_xlsx(list(sheet = p_net_pathway), paste0("../result/MPXV_VTG&IN&Wang_DEP_reactome_pathway_", toupper(reg), ".xlsx"))
}

# Plot (2022-12-16)
pathway_ca <- read.table('../data/Reactome_pathway/reactome_pathway_category.txt',
                         header = TRUE,
                         sep = '\t',
                         quote = '')
pathway_ca <- pathway_ca %>%
  rename(Description = pathway, Pathway_category = category)

p_pathway <- read_xlsx('../result/MPXV_VTG&IN&Wang_DEP_reactome_pathway_ALL.xlsx')
p_pathway <- p_pathway %>%
  select(-geneID) %>%
  left_join(pathway_ca, by = 'Description')
write_xlsx(list(sheet = p_pathway), '../result/MPXV_VTG&IN&Wang_DEP_reactome_pathway_ALL2.xlsx')

p_pathway2 <- p_pathway %>%
  group_by(Category) %>%
  top_n(n=-20, wt=p.adjust)
p_pathway2 <- rbind(p_pathway2, p_pathway[(p_pathway$Description == "Innate Immune System" | 
                                             p_pathway$Description == "Platelet activation, signaling and aggregation") &
                                            p_pathway$Category == "VTG",])
  
p_pathway2$Category <- factor(p_pathway2$Category, levels = c('VTG', 'IN', 'DEP'))

p_pathway2 <- p_pathway2 %>% arrange(Pathway_category, p.adjust)

p_pathway3 <- p_pathway2 %>%
  arrange(Pathway_category, p.adjust)

p_pathway_plot <- ggplot(p_pathway2, aes(Category, Description, color = Category)) +
  geom_point(aes(size = -log10(p.adjust))) +
  scale_y_discrete(limits = rev(unique(p_pathway3$Description))) +
  scale_x_discrete(limits = c('VTG', 'IN', 'DEP')) +
  scale_color_manual(values = c('#E74C3C', '#18BC9B', '#F1C50F')) +
  labs(x= '', y = '') +
  mytheme

ggsave(p_pathway_plot,
       filename = "../result/Fig_S4B.pdf",
       width = 6.5,
       height = 6.5)

vtp_pathway_ca_num <- p_pathway %>%
  filter(Category == 'VTG') %>%
  mutate_all(~replace(., is.na(.), 'Other')) %>%
  group_by(Pathway_category) %>%
  summarise(Number = n()) %>%
  arrange(desc(Number))

p_vtp_ca_num <- ggplot(vtp_pathway_ca_num, aes(Number, Pathway_category)) +
  geom_bar(stat = 'identity', fill = '#2880B9') +
  scale_y_discrete(limit = rev(vtp_pathway_ca_num$Pathway_category)) +
  scale_x_continuous(expand = c(0.03, 0.03)) +
  geom_text(aes(Number, label = Number)) +
  labs(y = 'Pathway Category') +
  mytheme
  
ggsave(p_vtp_ca_num,
       filename = "../result/Fig_S4C.pdf",
       width = 3,
       height = 5)
  
  