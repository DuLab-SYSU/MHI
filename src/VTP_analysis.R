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

library(gtable)
library(grid)

ggplot2.two_y_axis <- function(g1, g2) {
  g1 <- ggplotGrob(g1)
  g2 <- ggplotGrob(g2)
  
  pp <- c(subset(g1$layout, name == 'panel', se = t:r))
  
  g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == 'panel')]], 
                        pp$t, pp$l, pp$b, pp$l)
  
  hinvert_title_grob <- function(grob){
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  index <- which(g2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- g2$grobs[[index]]        # Extract that grob
  ylab <- hinvert_title_grob(ylab)     # Swap margins and fix justifications
  
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  index <- which(g2$layout$name == 'axis-l')  # Which grob
  yaxis <- g2$grobs[[index]]          # Extract the grob
  
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  yaxis$children[[2]] <- ticks
  
  g1 <- gtable_add_cols(g1, g2$widths[g2$layout[index, ]$l], pp$r)
  g1 <- gtable_add_grob(g1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(g1)
}

setwd('~/Desktop/Poxviridae-human/src/')
viruses <- c("Cowpox_virus", "Ectromelia_virus-ERPV", "Horsepox_virus", 
             "Molluscum_contagiosum_virus-subtype-1", "Monkeypox_virus", "Monkeypox_virus-Zaire-96-I-16",
             "Orf_virus", "Vaccinia_virus", "Vaccinia_virus-Ankara", 
             "Variola_virus", "Yaba-like_disease_virus", "Yaba_monkey_tumor_virus")
abbre <- c("CPXV", "ECTV ERPV", "HPXV", "MCV1", "MPXV", "MPXV Zaire-96-I-16",
           "ORFV", "VACV", "VACV Ankara", "VARV", "YLDV", "YMTV")

virus_abbre <- data.frame(Virus = viruses,
                          Abbreviation = abbre)

###########################################
### Network topological properties of VTGs 
### (2022-09-06)
node_centrality <- read.table('../intermediate/node_centrality.txt', 
                              header = TRUE, sep = '\t')


vtg_property <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  vtg <- read.table(paste0("../data/p-hipster/", v, ".txt"),
                    header = TRUE, sep = '\t')
  vtg <- dplyr::select(vtg, Human.Protein) %>%
    distinct() %>%
    dplyr::rename(Gene = Human.Protein)
  nodes <- unique(vtg$Gene)
  
  node_prop <- node_centrality %>%
    filter(Gene %in% nodes) %>%
    mutate(Virus = v) %>%
    group_by(Virus) %>%
    dplyr::summarise(Degree = mean(Degree),
                     Betweenness_centrality = mean(Betweenness_centrality),
                     Eigenvector_centrality = mean(Eigenvector_centrality),
                     Clustering_coefficient = mean(Clustering_coefficient),
                     Assortativity = mean(Assortativity),
                     Closeness_centrality = mean(Closeness_centrality))
  
  vtg_property <- rbind(vtg_property, node_prop)
}

randNones <- function(virus){
  
  vtg <- read.table(paste0("../data/p-hipster/", virus, ".txt"),
                           header = TRUE, sep = '\t')
  vtg <- dplyr::select(vtg, Human.Protein) %>%
    distinct() %>%
    dplyr::rename(Gene = Human.Protein)
  nodes <- unique(vtg$Gene)
  
  random_node <- data.frame(Number = 1:1000, Virus = virus)
  for (attr in colnames(node_centrality)[2:7]) {
    var <- assign(attr, c())
    
    for (i in 1:1000) {
      set.seed(i)
      node_samp <- node_centrality[sample(1:nrow(node_centrality), 
                                          size = length(nodes), replace = TRUE),]
      var <- c(var, mean(as.numeric(as.character(node_samp[attr][,1]))))
    }
    
    random_node[attr] <- var
  }
  
  return(random_node)
}

vtg_rand_property <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  
  random_node <- randNones(v)
  vtg_rand_property <- rbind(vtg_rand_property, random_node)
}

vtg_rand_property <- vtg_rand_property %>%
  inner_join(virus_abbre, by = 'Virus') %>%
  select(-Virus)

vtg_rand_property2 <- vtg_rand_property %>%
  group_by(Abbreviation) %>%
  dplyr::summarise(Degree_mean = mean(Degree),
                   Degree_sd = sd(Degree),
                   Betweenness_centrality_mean = mean(Betweenness_centrality),
                   Betweenness_centrality_sd = sd(Betweenness_centrality),
                   Eigenvector_centrality_mean = mean(Eigenvector_centrality),
                   Eigenvector_centrality_sd = sd(Eigenvector_centrality),
                   Clustering_coefficient_mean = mean(Clustering_coefficient),
                   Clustering_coefficient_sd = sd(Clustering_coefficient),
                   Assortativity_mean = mean(Assortativity),
                   Assortativity_sd = sd(Assortativity),
                   Closeness_centrality_mean = mean(Closeness_centrality),
                   Closeness_centrality_sd = sd(Closeness_centrality))

vtg_property <- inner_join(vtg_property, virus_abbre, by = 'Virus') %>% 
  select(-Virus) %>%
  left_join(vtg_rand_property2, by = "Abbreviation")
vtg_property2 <- vtg_property %>%
  mutate(Degree_Z = round((Degree - Degree_mean)/Degree_sd, 2),
         Degree_P = round(2*pnorm(-abs(Degree_Z)), 2),
         Betweenness_centrality_Z = round((Betweenness_centrality - 
                                       Betweenness_centrality_mean)/Betweenness_centrality_sd, 2),
         Betweenness_centrality_P = round(2*pnorm(-abs(Betweenness_centrality_Z)), 2),
         Eigenvector_centrality_Z = round((Eigenvector_centrality - 
                                       Eigenvector_centrality_mean)/Eigenvector_centrality_sd, 2),
         Eigenvector_centrality_P = round(2*pnorm(-abs(Eigenvector_centrality_Z)), 2),
         Clustering_coefficient_Z = round((Clustering_coefficient - 
                                       Clustering_coefficient_mean)/Clustering_coefficient_sd, 2),
         Clustering_coefficient_P = round(2*pnorm(-abs(Clustering_coefficient_Z)), 2),
         Assortativity_Z = round((Assortativity - Assortativity_mean)/Assortativity_sd, 2),
         Assortativity_P = round(2*pnorm(-abs(Assortativity_Z)), 2),
         Closeness_centrality_Z = round((Closeness_centrality - 
                                     Closeness_centrality_mean)/Closeness_centrality_sd, 2),
         Closeness_centrality_P = round(2*pnorm(-abs(Closeness_centrality_Z)), 2))

p_degree <- ggplot(vtg_rand_property, aes(Abbreviation, Degree)) +
  geom_boxplot(outlier.shape = 1) +
  geom_point(data = vtg_property2, aes(Abbreviation, Degree), size = 2, 
             color = "#E74C3C") +
  geom_text(data = vtg_property2, aes(x = Abbreviation, y = Degree + 0.5, label = Degree_Z)) +
  labs(x="",y="Degree") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

p_betweenness <- ggplot(vtg_rand_property, aes(Abbreviation, Betweenness_centrality)) +
  geom_boxplot(outlier.shape = 1) +
  geom_point(data = vtg_property2, aes(Abbreviation, Betweenness_centrality), size = 2, 
             color = "#E74C3C") +
  geom_text(data = vtg_property2, aes(x = Abbreviation, y = Betweenness_centrality + 0.00001, 
                                      label = Betweenness_centrality_Z)) +
  labs(x="",y="Betweenness centrality") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))

p_property <- ggpubr::ggarrange(p_degree, p_betweenness,
                                nrow=2, labels = c("A", "B"),
                                font.label = list(size = 11))
ggsave(p_property,
       filename = "../result/Fig_2AB.pdf",
       width = 9,
       height = 12,
       units = c("cm"))

#########################################################
### The viral module (the large connected component, LCC) 
### LCC size (2022-07-06)
lccPlot <- function(virus){
  expected_lcc_size <- read.table(paste0('../result/viral_module/',
                                         virus, '_expected_lcc_size.txt'),
                                  header = FALSE, sep = '\t')
  colnames(expected_lcc_size) <- c('LCC_size')
  lcc_statistic <-read.table(paste0('../result/viral_module/',
                                    virus, '_lcc_statistic.txt'),
                             header = TRUE, sep = '\t')
  
  p <- ggplot(expected_lcc_size, aes(LCC_size)) +
    geom_histogram(bins = 20, aes(y=..density..), fill = '#2880B9',
                   color = "white") +
    stat_function(aes(LCC_size),fun=dnorm, colour='#999999', 
                  xlim = c(min(expected_lcc_size$LCC_size), 
                           max(expected_lcc_size$LCC_size)),
                  args = list(mean = mean(expected_lcc_size$LCC_size), 
                              sd = sd(expected_lcc_size$LCC_size))) +
    geom_segment(aes(x = lcc_statistic$size, y=0.01,
                     xend = lcc_statistic$size, yend=0),
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "#E74C3C") +
    annotate("text", x = lcc_statistic$size-10, 
             y = 0.015, color = "#E74C3C",
             label = paste0(virus, " LCC\n", 
                            "Observed: ", lcc_statistic$size, 
                            "\nZ-Score: ", lcc_statistic$z.score), size = 2) +
    scale_y_continuous(expand = c(0.0005,0.0005)) +
    labs(x = 'LCC', y = 'Density') +
    mytheme
  
  return(p)
}

for (f in list.files(path = "../data/p-hipster/")){
  virus <- sub(".txt", "", f)
  
  p <- lccPlot(virus)
  # Figure S1
  ggsave(p,
         filename = paste0('../intermediate/viral_module/',
                           virus, '_LCC.pdf'),
         width = 5,
         height = 5,
         units = c("cm"))
}

### % targets formed LCC (2022-09-05)
library(igraph)
ppi <- read.table('../result/integrated_ppi_largest_component.txt',
                  header = FALSE,
                  sep = '\t',
                  quote = '')
ppin <- graph_from_data_frame(ppi, directed = FALSE)
nodes <- get.vertex.attribute(ppin)
nodes <- nodes$name

ptg <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  vtg <- read.table(paste0("../data/p-hipster/", v, ".txt"),
                    header = TRUE, sep = '\t')
  vtg <- dplyr::select(vtg, Virus.Protein, Human.Protein) %>%
    distinct() %>%
    dplyr::rename(Gene = Human.Protein, Viral_protrin = Virus.Protein) %>%
    mutate(Virus = v) %>%
    filter(Gene %in% nodes)
  
  ptg <- rbind(ptg, vtg)
}

lcc_perc <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  lcc_stat <- read.table(paste0("../intermediate/viral_module/", v, "_lcc_statistic.txt"),
                    header = TRUE, sep = '\t')
  lcc_stat <- lcc_stat %>% 
    mutate(Virus = v)
  
  lcc_perc <- rbind(lcc_perc, lcc_stat)
}

lcc_perc <- ptg %>% 
  dplyr::select(Gene, Virus) %>% 
  distinct() %>% 
  group_by(Virus) %>% 
  dplyr::summarise(Target_num = n()) %>%
  inner_join(lcc_perc, by = 'Virus') %>%
  mutate(LCC_percent = round(size/Target_num, 4)*100)

lcc_perc <- inner_join(lcc_perc, virus_abbre, by = 'Virus')

p_lcc_size <- ggplot(lcc_perc, aes(Abbreviation, size)) +
  geom_bar(stat = 'identity', fill = '#2880B9', alpha = .75) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x="",y="LCC size") +
  mytheme +
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

p_lcc_perc <- ggplot(lcc_perc, aes(Abbreviation, LCC_percent)) +
  geom_line(group = 1, color = '#E74C3C', size = 1.2) +
  geom_point(color = '#E74C3C', size = 2) +
  # geom_text(data = lcc_perc, 
  #           aes(x = Abbreviation, y = size - 100, label = LCC_percent),
  #           size = 2.5, color = "white") +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(90, 96)) +
  labs(x="",y="LCC proportion") +
  mytheme +
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1))

# Figure 2C
ggplot2.two_y_axis(p_lcc_size, p_lcc_perc)

# ggsave(p_lcc_perc,
#        filename = "../result/viral_LCC2.pdf",
#        width = 9,
#        height = 8,
#        units = c("cm"))

#################################
### Pathway module separation 
### (2022-09-05)
pathway_separation_z <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  ps_z <- read.table(paste0("../intermediate/pathway_separation/", v, "_separation_zscore.txt"),
                   header = TRUE, sep = '\t', quote = '')
  ps_z <- ps_z %>% mutate(Virus = v)
  
  pathway_separation_z <- rbind(pathway_separation_z, ps_z)
}
pathways <- unique(pathway_separation_z$pathway)

pathway_separation <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  ps <- read.table(paste0("../intermediate/pathway_separation/", v, "_separation.txt"),
                     header = TRUE, sep = '\t', quote = '')
  ps <- ps %>% 
    mutate(Virus = v) %>%
    filter(pathway %in% pathways)
  
  pathway_separation <- rbind(pathway_separation, ps)
}

pathway_separation2 <- pathway_separation %>%
  inner_join(virus_abbre, by = 'Virus') %>%
  select(pathway, separation, Abbreviation) %>%
  pivot_wider(names_from = Abbreviation, values_from = separation) %>%
  arrange(pathway)

pathway_separation_z2 <- pathway_separation_z %>%
  inner_join(virus_abbre, by = 'Virus') %>%
  dplyr::select(Abbreviation, pathway, z.score) %>%
  pivot_wider(names_from = Abbreviation, values_from = z.score) %>%
  arrange(pathway)
pathway_separation_z2[is.na(pathway_separation_z2)] <-  0

library(ComplexHeatmap)
library(circlize)

pdf('../result/Fig_2D.pdf',
    width = 9,
    height = 10)
Heatmap(as.matrix(pathway_separation2[,2:13]),
        row_labels = pathway_separation2$pathway,
        cluster_columns = FALSE,
        col = colorRamp2(c(min(as.matrix(pathway_separation2[,2:13])), 0, 
                           max(as.matrix(pathway_separation2[,2:13]))), 
                         c("#3398DA", "white", "#E74C3C")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ifelse(as.matrix(pathway_separation_z2[,2:13])[i, j] < -1.6, "*", ""),
                    x, y, gp = gpar(fontsize = 10))
        })
dev.off()


#########################
### Disease Comorbidities
### (2022-09-09)
disease_separation <- data.frame()
for (f in list.files(path = "../data/p-hipster/")){
  v <- sub(".txt", "", f)
  ds <- read.table(paste0("../intermediate/disease_separation/", v, "_separation.txt"),
                          header = TRUE, sep = '\t', quote = '')
  ds <- ds %>% mutate(Virus = v)
  
  disease_separation <- rbind(disease_separation, ds)
}

### Pearson correlation coefficients of the separation for each viral pair (2022-09-09)
disease_separation_wide <- inner_join(disease_separation, virus_abbre, by = 'Virus') %>%
  select(disease, separation, Abbreviation) %>%
  pivot_wider(names_from = Abbreviation, values_from = separation)

separation_pearson <- as.data.frame(t(combn(abbre,2)))
s_pear_cor <- function(separation_pearson, c1, c2){
  values <- cor.test(as.numeric(as.matrix(disease_separation_wide[separation_pearson[c1]])), 
                     as.numeric(as.matrix(disease_separation_wide[separation_pearson[c2]])), 
                     method="pearson")
  values$estimate
}


s_pear_p <- function(separation_pearson, c1, c2){
  values <- cor.test(as.numeric(as.matrix(disease_separation_wide[separation_pearson[c1]])), 
                     as.numeric(as.matrix(disease_separation_wide[separation_pearson[c2]])), 
                     method="pearson")
  values$p.value
}


colnames(separation_pearson) <- c("virus_A", "virus_B")
separation_pearson$correlation <- apply(separation_pearson[,1:2], 1, 
                                        s_pear_cor, c1 = 'virus_A', c2 = 'virus_B')
separation_pearson$p_value <- apply(separation_pearson[,1:2], 1, 
                                    s_pear_p, c1 = 'virus_A', c2 = 'virus_B')

separation_pearson2 <- separation_pearson %>%
  dplyr::select(virus_B, virus_A, correlation, p_value) %>%
  dplyr::rename(virus_A = virus_B, virus_B = virus_A) %>%
  rbind(separation_pearson) %>%
  dplyr::select(-p_value) %>%
  pivot_wider(names_from = virus_B, values_from = correlation) %>%
  arrange(virus_A)
separation_pearson2[is.na(separation_pearson2)] <- 1

pdf('../result/Fig_3B.pdf',
    width = 10,
    height = 7)
Heatmap(as.matrix(separation_pearson2[,2:13]),
        name = 'Pearson correlation coefficient',
        row_labels = separation_pearson2$virus_A,
        col = colorRamp2(c(min(as.matrix(separation_pearson2[,2:13])), 
                           max(as.matrix(separation_pearson2[,2:13]))), 
                         c("white", "#E74C3C")))
dev.off()

### Disease comorbidity measured by the network overlap between MPXV targets and 299 diseases (2022-09-12)
## disease category and gene number for each disease
disease_gene <- read.table('../data/disease_genes/Menche2015/disease_genes.txt',
                           header = TRUE, sep = '\t', quote = '') 
disease_gene_num <- filter(disease_gene, gene %in% nodes) %>%
  distinct() %>%
  group_by(disease) %>%
  summarise(number = n())
disease_category <- read.table('../data/disease_genes/Menche2015/disease_category.txt',
                               header = TRUE, sep = '\t', quote = '')
disease_category <- disease_category %>%
  inner_join(disease_gene_num, by = 'disease') %>%
  separate_rows(category, sep = ';')

## Network separation between MPXV targets and 299 diseases (2023-01-11)
colors <- c('#10A185', '#27AE61', '#7F8D8D', '#8E45AD', '#2D3D50',
            '#F39C12', '#D35301', '#BF3A2B', '#2880B9', '#BDC3C8')

MPXV_disease_separation <- inner_join(disease_separation, virus_abbre, by = 'Virus') %>%
  filter(Abbreviation == 'MPXV') %>%
  select(-Virus) %>%
  inner_join(disease_category, by = 'disease') %>%
  arrange(category)

to_add <- data.frame(matrix(NA, 5, ncol(MPXV_disease_separation)))
colnames(to_add) <- colnames(MPXV_disease_separation)
to_add2 <- data.frame(matrix(NA, 30, ncol(MPXV_disease_separation)))
colnames(to_add2) <- colnames(MPXV_disease_separation)


MPXV_disease_separation2 <- rbind(to_add, MPXV_disease_separation, to_add2) %>%
  mutate(id = seq(1, nrow(MPXV_disease_separation)+35, 1))
label_data <- MPXV_disease_separation2 %>% 
  drop_na() %>% 
  filter(separation < 0)
label_data2 <- MPXV_disease_separation2 %>% 
  drop_na() %>%
  group_by(category) %>% 
  summarise(mean = mean(id))

library(ggrepel)
p_sepa <- ggplot() +
  geom_segment(aes(x = 2, y = 0, xend = 372, yend = 0),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 2, y = 0.25, xend = 372, yend = 0.25),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 2, y = 0.5, xend = 372, yend = 0.5),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 2, y = 0.75, xend = 372, yend = 0.75),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 2, y = 1, xend = 372, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 2, y = 0, xend = 2, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 76, y = 0, xend = 76, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 150, y = 0, xend = 150, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 224, y = 0, xend = 224, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 298, y = 0, xend = 298, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_segment(aes(x = 372, y = 0, xend = 372, yend = 1),
               color = '#BDC3C8', alpha = 1, size = .3, inherit.aes = FALSE) +
  geom_point(MPXV_disease_separation2, mapping = aes(id, separation, color = category, size = number)) +
  geom_text_repel(label_data, mapping = aes(id, separation, label = disease), size = 3) +
  geom_text(label_data2, mapping = aes(mean, 1.1, label = category, color = category), size = 3) +
  annotate("text", x = max(MPXV_disease_separation2$id), y = c(0, 0.25, 0.50, 0.75, 1),
           label = c("0.00", "0.25", "0.50", "0.75", "1.00"), 
           color = "grey", size = 3, angle = 0, fontface = "bold", hjust = 1) +
  scale_color_manual(values = colors) +
  labs(x = '', y = 'Separation') +
  ylim(-1.2, 1.2) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  coord_polar()
  
ggsave(p_sepa,
       filename = "../result/Fig_3A.pdf",
       width = 18,
       height = 18,
       units = c("cm"))


##################################################
### MPXV targets, IN and DEPs after MPXV infection

## MPXV targets, IN and DEPs after MPXV infection
vtgs <- read.table('../data/p-hipster/Monkeypox_virus.txt',
                   header = TRUE,
                   sep = '\t',
                   quote = '')
vtgs <- unique(vtgs$Human.Protein)

# Plasma specific network 
plasma_ppi <- read.table('../result/plasma_specific_network.txt',
                  header = FALSE,
                  sep = '\t',
                  quote = '')
plasma_net <- graph_from_data_frame(plasma_ppi, directed = FALSE)
pnodes <- get.vertex.attribute(plasma_net)
pnodes <- pnodes$name

p_vtgs <- vtgs[vtgs %in% pnodes]
## Gene overlap (2022-10-28)
deps <- read.table('../data/proteome/MPXV_plasma_DEP.txt',
                   header = FALSE,
                   sep = '\t',
                   quote = '')

deps <- filter(deps, V1 %in% pnodes)
all_deps <- unique(deps$V1)
up_deps <- deps[deps$V2=='up',1]
down_deps <- deps[deps$V2=='down',1]

ins <- read.table('../intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/intermediate_nodes_all.txt',
                  header = FALSE,
                  sep = '\t')
ins <- ins$V1

pdf(file = "../result/Fig_S4A.pdf")
tid_plot <- Venn(list("VTGs" = p_vtgs,
                      "INs" = ins,
                      "DEPs" = all_deps))
plot(tid_plot, doWeight = T)
dev.off()










