Biomart_geneid.txt # ID mapping from biomart Ensembl genome browser 104 (2021-11-16)

### poxviruses-human PPIs
p-hipster/ # phipster.org (2022-07-23)

### human PPIs
cd PPIs_from_databases/
BIOGRID-ORGANISM-Homo_sapiens-4.4.203.mitab.txt.gz # https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.203/ (2021-11-09)
Hsapi20170205.txt.gz # https://dip.doe-mbi.ucla.edu/dip/ (2021-11-09)
MINT_human.gz # http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/current/search/query/species:human (2021-11-09)
intact.txt.gz # ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip (2021-11-09)
matrixdb_FULL.tab.gz # http://matrixdb.univ-lyon1.fr (2022-03-09)


### PSI-MI
cd ../../src/PSI-MI/
mi.owl # relationship between MI, https://www.ebi.ac.uk/ols/ontologies/MI (2021-11-10)

### Reactome pathways
cd ../../data/Reactome_pathway
c2.cp.reactome.v7.5.1.symbols.gmt # http://www.gsea-msigdb.org/gsea/downloads.jsp (2022-07-12)
reactome_pathways_from_MsigDB.txt # (src/reactome_pathway.ipynb) (2022-07-13)
UniProt2Reactome_All_Levels.txt # All levels of the pathway hierarchy (2022-12-15) 
reactome_pathway_category.txt # (src/reactome_pathway.ipynb) (2022-12-16)

### disease genes
cd ../disease_genes/Menche2015
disease_genes.txt # disease-associated genes from Menche et al., Science 2015
disease_category.txt # disease categories, https://www.ncbi.nlm.nih.gov/mesh/ (2022-09-07)

### Proteomic data
cd ../../proteome
proteome/human_plasma_proteome.xlsx # Human Plasma Proteome Project Data from Deutsch et al. 2021 (www.peptideatlas.org/hupo/hppp/) (2022-10-25)
MPXV_plasma_DEP.txt # The differential expressed proteins after monkeypox infection from Wang et al. 2022 (2022-10-26)

### Drug-target interactions
cd ../DrugBank
full_database.xml # DrugBank v5.1.9 https://go.drugbank.com/releases/latest (2022-06-26)
../../src/parse_drugbank.R # (2022-07-04)

### side effect
cd ../SIDER # Side Effect Resource (SIDER 4.1) sideeffects.embl.de/download/ (2022-11-15)
../../src/side_effect_proteins.R


### REMAP/
## protein sequences
uniprot-download_1to18101.fasta # sequences for 18101 proteins extracted from UniProt (2022-10-21)
uniprot_human_proteins.fas #  rename the header name (protein_index.ipynb) (2022-10-21)
cat uniprot_human_proteins.fas | sed '/^>/ s/$/ /' | tr -s '\s' | tr -d "\n" | sed 's/>/\n>/g' | sed '1 d' | awk '{print $1, length($2), $2}' | sort -k1,1V -k2nr | grep -v "^>no|" | awk 'L!=$1 {print $1"\n"$3}{L=$1}' > uniprot_human_proteins_noRedu.fas # no redundant protein sequences (2022-10-21)

MPXV_protein_seq.faa # all MPXV protein sequences, https://download.cncb.ac.cn/Genome/Viruses/Poxviridae/Monkeypox_virus/ (2022-11-28)
cat MPXV_protein_seq.faa | sed '/^>/ s/$/ /' | tr -s '\s' | tr -d '\n' | sed 's/>/\n>/g' | sed '1 d' | grep "\[Monkeypox virus\]" | sed 's/ |.*complete//' | sed 's/>/>MPXV|MPXV-/' | awk '{print $1, length($2), $2}' | sort -k1,1V -k2nr | awk 'L!=$1 {print $1"\n"$3}{L=$1}' > MPXV_protein_seq2.faa
cd-hit -i MPXV_protein_seq2.faa -o MPXV_protein_noRedu.fas -c 0.95 -T 24 # identification of no redundant protein sequences by cd-hit (2022-11-28)

cat uniprot_human_proteins_noRedu.fas MPXV_protein_noRedu.fas > all_proteins.fas # all human and MPXV protein sequences (2022-11-28)

