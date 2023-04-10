library(tidyverse)
library(RColorBrewer)
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
library(igraph)
ppi <- read.table('../result/integrated_ppi_largest_component.txt',
                  header = FALSE,
                  sep = '\t',
                  quote = '')
ppin <- graph_from_data_frame(ppi, directed = FALSE)
nodes <- get.vertex.attribute(ppin)
nodes <- nodes$name

dti <- read.table('../data/REMAP/drug_targets_approved.tsv',
                  header = TRUE,
                  sep = '\t',
                  quote = '') %>%
  rename(protein = gene)

raw_pred <- read.table('../data/REMAP/REMAP_pred.tsv',
                       header = TRUE,
                       sep = '\t',
                       quote = '')
raw_pred <- dti %>%
  mutate(label = 'active') %>%
  right_join(raw_pred, by = c('drug', 'protein')) %>%
  mutate_all(~replace(., is.na(.), 'inactive')) %>%
  mutate(score = as.numeric(score)) %>%
  filter(score <= 1.1)

obs_dti <- raw_pred %>%
  filter(label == 'active')

p_dti_score <- ggplot(obs_dti, aes(score)) +
  geom_density() +
  labs(x = 'Raw prediction score', y = 'Density') +
  mytheme

ggsave(p_dti_score,
       filename = "../result/Fig_S2A.pdf",
       width = 8,
       height = 6,
       units = c("cm"))

### REMAP score distributions (2022-11-24)
obs_dti_num <- obs_dti %>%
  select(score) %>%
  mutate(interval = cut_width(score, width = 0.05, center = 0.025)) %>%
  count(interval) %>%
  rename(obs_num = n)

pred_dti_num <- raw_pred %>%
  select(score) %>%
  mutate(interval = cut_width(score, width = 0.05, center = 0.025)) %>%
  count(interval) %>%
  rename(all_num = n) %>%
  left_join(obs_dti_num, by = 'interval') %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(freq = round(obs_num/all_num, 4)*100) %>%
  mutate(center = seq(0.025, 1.075, 0.05))


### Performance of REMAP (2022-11-25)
library(foreach)
library(doParallel)

#st = Sys.time()
cl <- makeCluster(25)
registerDoParallel(cl)

raw_pred_top1perc <- raw_pred %>%
  filter(score > quantile(raw_pred$score, 0.99)) %>%
  arrange(desc(score))

TPRs <- c()
tprFun <- function(x){
  library(dplyr)
  rank_score <- raw_pred_top1perc[1:x,] %>%
    filter(label == 'active')
  TPR <- nrow(rank_score)/nrow(obs_dti)
  TPR
}

TPRs <- foreach(i = 1:nrow(raw_pred_top1perc),
                .combine = c) %dopar%
  tprFun(i)

stopCluster(cl)

rank_TPR <- data.frame(cutoff_rank = 1:nrow(raw_pred_top1perc),
                       TPR = TPRs)
write.table(rank_TPR,
            '../result/REMAP/REMAP_score_cutoffRank_TPR.txt',
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

rank_TPR <- read.table('../result/REMAP/REMAP_score_cutoffRank_TPR.txt',
                       header = TRUE,
                       sep = '\t')

### Cutoff at 80% of the observed DTIs were included in predicted DTIs (2022-12-05)
raw_pred_80perc <- raw_pred %>%
  filter(score >= quantile(obs_dti$score, 0.2))

write.table(raw_pred_80perc,
            '../result/REMAP/REMAP_pred30256.txt',
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

p_rank_TPR <- ggplot(rank_TPR, aes(cutoff_rank, TPR)) +
  geom_line() +
  geom_point(data = rank_TPR[rank_TPR$cutoff_rank==30256,], aes(cutoff_rank, TPR), color = "#E74C3C") +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = 'CutoffRank', y = 'True Positive Rate') +
  annotate("text", x = 65000, y = 0.75, color = "#E74C3C",
           label = "Cutoff: 30,256\nTPR: 0.8") +
  mytheme

ggsave(p_rank_TPR,
       filename = "../result/Fig_S2B.pdf",
       width = 8,
       height = 6,
       units = c("cm"))



