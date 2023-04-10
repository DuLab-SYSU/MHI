library(dbparser)
library(tidyverse)
library(ggplot2)
library(XML)

setwd('~/Documents/Computational_Biology/Poxviridae-human/src/')
### (2022-07-04)
## parse data from XML and save it to memory
read_drugbank_xml_db("~/Documents/Computational_Biology/Poxviridae-human/data/DrugBank/full_database.xml")

## load drugs general information
drugs <- drug_general_information()
drugs <- dplyr::select(drugs, primary_key, type, name, state) %>%
  dplyr::rename(ID = primary_key, Name = name, State = state)

## load drug groups data
drug_groups <- drug_groups()
drugs <- inner_join(drugs, drug_groups, by = c("ID" = "drugbank-id"))

## load drug ATC code data
drug_atc_codes <- drug_atc_codes()
#writexl::write_xlsx(list(drug_atc_code = drug_atc_codes),
#           '../data/DrugBank/drug_ATC_codes.xlsx')

drug_atc_codes2 <- drug_atc_codes %>%
  dplyr::select(atc_code, `drugbank-id`) %>%
  dplyr::rename(ATC_code = atc_code, ID = `drugbank-id`)
drug_atc_codes2 <- aggregate(drug_atc_codes2[1], drug_atc_codes2[2],
                             FUN = function(X) paste(unique(X), collapse=";"))

## load drug targets data
drug_targets <- targets_polypeptides()
drug_targets <- drug_targets %>% 
  filter(organism == "Humans") %>%
  dplyr::select(name, gene_name, cellular_location, parent_id)
drug_targets2 <- targets()
drug_targets <- dplyr::select(drug_targets2, id, parent_key) %>%
  inner_join(drug_targets, by = c("id" = "parent_id")) %>%
  filter(parent_key %in% drugs$ID) %>%
  mutate(Type = "polypeptide")

## load enzymes data
drug_enzymes <- enzymes_polypeptides()
drug_enzymes <- drug_enzymes %>% 
  filter(organism == "Humans") %>%
  dplyr::select(name, gene_name, cellular_location, parent_id)
drug_enzymes2 <- enzymes()
drug_enzymes <- dplyr::select(drug_enzymes2, id, parent_key) %>%
  inner_join(drug_enzymes, by = c("id" = "parent_id")) %>%
  filter(parent_key %in% drugs$ID) %>%
  mutate(Type = "enzyme")

## load carriers data
drug_carriers <- carriers_polypeptides()
drug_carriers <- drug_carriers %>% 
  filter(organism == "Humans") %>%
  dplyr::select(name, gene_name, cellular_location, parent_id)
drug_carriers2 <- carriers()
drug_carriers <- dplyr::select(drug_carriers2, id, parent_key) %>%
  inner_join(drug_carriers, by = c("id" = "parent_id")) %>%
  filter(parent_key %in% drugs$ID) %>%
  mutate(Type = "carrier")

## load transporters data
drug_transporters <- transporters_polypeptides()
drug_transporters <- drug_transporters %>% 
  filter(organism == "Humans") %>%
  dplyr::select(name, gene_name, cellular_location, parent_id)
drug_transporters2 <- transporters()
drug_transporters <- dplyr::select(drug_transporters2, id, parent_key) %>%
  inner_join(drug_transporters, by = c("id" = "parent_id")) %>%
  filter(parent_key %in% drugs$ID) %>%
  mutate(Type = "transporter")


### drugs and targets
final_drug_target <- rbind(drug_targets,
                           drug_enzymes,
                           drug_carriers,
                           drug_transporters) %>%
  dplyr::rename(DB_id = id, ID = parent_key, Gene_target = gene_name) %>%
  inner_join(drugs, by = "ID") %>%
  left_join(drug_atc_codes2, by = "ID") %>%
  distinct()

final_drug_target <- final_drug_target %>%
  dplyr::select(type, ID, Name, State, group, Gene_target,
                DB_id, name, Type, ATC_code) %>%
  filter(Gene_target != "")
write.table(final_drug_target, 
            "~/Documents/Computational_Biology/Poxviridae-human/data/DrugBank/drug_targets_full.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

approved_drug_target <- final_drug_target %>%
  filter(group == "approved")
write.table(approved_drug_target, 
            "~/Documents/Computational_Biology/Poxviridae-human/data/DrugBank/drug_targets_approved.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

polypeptide_drug_target <- approved_drug_target %>%
  filter(Type == "polypeptide")
write.table(polypeptide_drug_target, 
            "~/Documents/Computational_Biology/Poxviridae-human/data/DrugBank/drug_targets_polypeptide.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

### Number of drug and target (2022-08-02)
setwd('~/Documents/Computational_Biology/Poxviridae-human/data/DrugBank/')
ppi <- read.table('../integrated_ppi_largest_component.txt',
                  header = FALSE,
                  sep = '\t')
nodes <- unique(c(as.character(ppi$V1), as.character(ppi$V2)))
polypeptide_drug_target <- read.table('drug_targets_polypeptide.txt',
                                      header = TRUE,
                                      sep = '\t',
                                      quote = '')
drug_target_in_interactome <- polypeptide_drug_target %>%
  filter(Gene_target %in% nodes)
  
### Drug ATC codes (2022-09-10)
drug_targets <- read.table('../data/DrugBank/drug_targets_full.txt',
                           header = TRUE, sep = '\t', quote = '')
drugs <- select(drug_targets, ID, Name, ATC_code) %>%
  rename(atc_code = ATC_code, `drugbank-id` = ID, name = Name) %>%
  distinct() %>%
  separate_rows(atc_code, sep = ';')

drug_atc_codes <- drugs %>%
  left_join(drug_atc_codes, by = c('atc_code', 'drugbank-id'))

writexl::write_xlsx(list(drug_atc_code = drug_atc_codes),
                    '../data/DrugBank/drug_ATC_codes.xlsx')

## ATC code level4 (2022-09-10)
drug_atc_codes <- read_xlsx('../data/DrugBank/drug_ATC_codes.xlsx')
code_level4 <- drug_atc_codes %>%
  select(name, level_4, code_4) %>%
  distinct() %>%
  mutate_at(c(2:3), ~replace(., is.na(.), 'Others')) %>%
  mutate(level_4 = str_to_sentence(level_4))

code_level4_2 <- aggregate(code_level4[2], code_level4[1],
                           FUN = function(X) paste(unique(X), collapse=";"))
code_level4_2$level_4_category <- apply(code_level4_2[,2,drop = FALSE], 1, 
                                        function(x) ifelse(length(str_split(x, ';')[[1]]) > 1, 'Various', x))

writexl::write_xlsx(list(code_level4 = code_level4_2),
                    '../data/DrugBank/drug_ATC_codes_level4.xlsx')

### SMILE (2022-10-24)
approved_drug_target <- read.table('../data/DrugBank/drug_targets_approved.txt',
                                   header = TRUE, sep = '\t', quote = '')

library(readxl)
drug_smi <- read_xlsx('../data/DrugBank/structure links.xlsx')
drug_smi2 <- drug_smi %>%
  rename(ID = `DrugBank ID`) %>%
  select(ID, Name, InChIKey, SMILES) %>%
  filter(Name %in% unique(approved_drug_target$Name)) %>%
  drop_na()
write.table(drug_smi2, '../data/DrugBank/drug_structure_approved.tsv',
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)







