#!/usr/bin/env bash

cd ..
data=$PWD/data
intermediate=$PWD/intermediate
result=$PWD/result
src=$PWD/src

###################
### node properties
python $src/node_property.py \
    -i $data/integrated_ppi_largest_component.txt \
    -o $intermediate/node_property.txt


#################
###  viral module
mkdir -p $intermediate/viral_module

num_cores=10

parallel -u -j $num_cores --bar \
    python $src/viral_module.py \
        -v $data/p-hipster/{}.txt \
        -i $result/integrated_ppi_largest_component.txt \
        -o1 $intermediate/viral_module/{}_lcc_statistic.txt \
        -o2 $intermediate/viral_module/{}_expected_lcc_size.txt \
        ::: `ls $data/p-hipster/ | sed 's/.txt//'`
        
######################
### pathway separation
mkdir -p $intermediate/pathway_separation

for virus in `ls $data/p-hipster/ | sed 's/.txt//'`
do
    echo $virus

    python $src/pathway_separation.py \
        -v $data/p-hipster/"$virus".txt \
        -p $data/Reactome_pathway/reactome_pathways_from_MsigDB.txt \
        -i $result/integrated_ppi_largest_component.txt \
        -o $intermediate/pathway_separation/"$virus"_separation.txt \
        -o2 $intermediate/pathway_separation/"$virus"_separation_zscore.txt
done
        

######################
### disease separation
mkdir -p $result/disease_separation

for virus in `ls $data/p-hipster/ | sed 's/.txt//'`
do
    echo $virus

    python $src/disease_separation.py \
        -v $data/p-hipster/"$virus".txt \
        -d $data/disease_genes/Menche2015/disease_genes.txt \
        -i $result/integrated_ppi_largest_component.txt \
        -o $intermediate/disease_separation/"$virus"_comorbidity.txt
done


######################
### intermediate nodes
mkdir -p $intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022
for reg in 'all' 'up' 'down'
do
    echo $reg   	
    /share/home/dxj_tangkang/software/Python-3.9.12/python $src/intermediate_nodes.py \
        -v $data/p-hipster/Monkeypox_virus.txt \
        -d $data/proteome/MPXV_plasma_DEP.txt \
        -r $reg \
        -i $result/plasma_specific_network.txt \
        -o $intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/intermediate_nodes_"$reg".txt \
        -o2 $intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/significant_paths_"$reg".txt
done

### Reactome pathway enrichment for 3 protein sets
Rscript --vanilla $src/pathway_enrichment.R

#################################
### analysis visualization of VTPs

Rscript --vanilla $src/VTP_analysis.R


#########################################
### drug proximity based on DrugBank DTIs

### potential drugs
for virus in `ls $data/p-hipster/ | sed 's/.txt//'`
do
    echo $virus
    mkdir -p $intermediate/DrugBank/drug_proximity/$virus

    python3 $src/drug_proximity.py \
        -v $data/p-hipster/"$virus".txt \
        -d $data/DrugBank/drug_targets_approved.txt \
        -i $result/integrated_ppi_largest_component.txt \
        -o $intermediate/DrugBank/drug_proximity/$virus
done

# statistics of drug candidates
Rscript --vanilla $src/DrugBank_drugs.R

### network-based drug and pathway proximity
mkdir -p $intermediate/DrugBank/drug_pathway_proximity

python3 $src/drug_pathway_proximity.py \
    -p $data/Reactome_pathway/reactome_pathways_from_MsigDB.txt \
    -d $data/DrugBank/drug_targets_approved.txt \
    -l $result/MPXV_277candidate_drugs.txt \
    -i $result/integrated_ppi_largest_component.txt \
    -o $intermediate/DrugBank/drug_pathway_proximity/MPXV_277candidate_drugs_pathway_proximity.txt
    
    
### network-based drug and side effect proximity
mkdir -p $intermediate/DrugBank/drug_sideEffect_proximity

python3 $src/drug_sideEffect_proximity.py \
    -p $data/SIDER/side_effect_proteins.txt \
    -d $data/DrugBank/drug_targets_approved.txt \
    -l $result/MPXV_277candidate_drugs.txt \
    -i $result/integrated_ppi_largest_component.txt \
    -o $intermediate/DrugBank/drug_sideEffect_proximity/MPXV_277candidate_drugs_sideEffect_proximity.txt


### network proximity between drug targets and 3 protein sets

## distance between drug targets and MPXV targets (plasma specific network)
mkdir -p $intermediate/DrugBank/drug_proximity_plasma/VTP

python $src/drug_proximity.py \
    -v $data/p-hipster/Monkeypox_virus.txt \
    -d $data/DrugBank/drug_targets_approved.txt \
    -i $result/plasma_specific_network.txt \
    -o $intermediate/DrugBank/drug_proximity_plasma/VTP


## distance between drug targets and DEPs
for reg in 'all' 'up' 'down'
do
    echo $reg
    mkdir -p $intermediate/DrugBank/drug_proximity_plasma/DEP/plasma_proteome_from_Wang2022/$reg
    
    python $src/drug_proximity_forDEPs.py \
        -p $data/proteome/MPXV_plasma_DEP.txt \
        -d $data/DrugBank/drug_targets_approved.txt \
        -r $reg \
        -i $result/plasma_specific_network.txt \
        -o $intermediate/DrugBank/drug_proximity_plasma/DEP/plasma_proteome_from_Wang2022/$reg
done

## distance between drug targets and INs
for reg in 'all' 'up' 'down'
do
    echo $reg
    mkdir -p $intermediate/DrugBank/drug_proximity_plasma/IN/$reg

    python $src/drug_proximity_forINs.py \
        -g $intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/intermediate_nodes_"$reg".txt \
        -d $data/DrugBank/drug_targets_approved.txt \
        -i $result/plasma_specific_network.txt \
        -o $intermediate/DrugBank/drug_proximity_plasma/IN/$reg
done


#########################################
### drug proximity based on REMAP DTIs

### systematical prediction of DTIs
cd $data
# Protein-protein similatity (2022-10-22)
python $src/REMAP/get_blast_sim.py \ 
    --input uniprot_human_proteins_noRedu.fas \
    --output human_prot_prot_sim.tsv

# Tanimoto 2D chemical similarity for drug pairs (2022-11-20)
python $src/REMAP/get_tani_sim.py \
    -s $data/REMAP/drug_structure_approved.tsv \
    -o $data/REMAP/drug_drug_sim.tsv

## REMAP prediction (2022-11-22)
python REMAP.py \
    --path $data \
    --R drug_targets_approved.tsv \
    --chemsim drug_drug_sim.tsv \
    --protsim human_prot_prot_sim.tsv \
    --seed 19110935 \
    --out_file $data/REMAP/REMAP_pred.tsv

### cutoff selection
Rscript --vanilla $src/REMAP/cutoff_selection.R

### potential drugs
for virus in `ls $data/p-hipster/ | sed 's/.txt//'`
do
    echo $virus
    mkdir -p $intermediate/REMAP/drug_proximity/$virus

    python $src/REMAP/drug_proximity.py \
        -v $data/p-hipster/"$virus".txt \
        -d $result/REMAP/REMAP_pred30256.txt \
        -i $result/integrated_ppi_largest_component.txt \
        -o $intermediate/REMAP/drug_proximity/$virus
done

# statistics of drug candidates
Rscript --vanilla $src/REMAP/REMAP_drugs.R

### network-based drug and side effect proximity
mkdir -p $intermediate/REMAP/drug_sideEffect_proximity

python $src/REMAP/drug_sideEffect_proximity.py \
    -p $data/SIDER/side_effect_proteins.txt \
    -d $result/REMAP/REMAP_pred30256.txt \
    -l $result/REMAP/candidate_drugs.txt \
    -i $result/integrated_ppi_largest_component.txt \
    -o $intermediate/REMAP/drug_sideEffect_proximity/

# comparison of side-effects of candidate drugs between DrugBank and REMAP DTIs
Rscript --vanilla $src/REMAP/REMAP_DrugBank_comparison.R


### network proximity between drug targets and 3 protein sets

## distance between drug targets and MPXV targets (plasma specific network)
mkdir -p $intermediate/REMAP/drug_proximity_plasma/VTP

python $src/REMAP/drug_proximity.py \
    -v $data/p-hipster/Monkeypox_virus.txt \
    -d $result/REMAP/REMAP_pred30256.txt \
    -i $result/plasma_specific_network.txt \
    -o $intermediate/REMAP/drug_proximity_plasma/VTP
    
## distance between drug targets and DEPs (plasma specific network)
mkdir -p $intermediate/REMAP/drug_proximity_plasma/DEP/
    
python $src/REMAP/drug_proximity_forDEPs.py \
    -p $data/proteome/MPXV_plasma_DEP.txt \
    -n $intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/intermediate_nodes_all.txt \
    -d $result/REMAP/REMAP_pred30256.txt \
    -i $result/plasma_specific_network.txt \
    -o $intermediate/REMAP/drug_proximity_plasma/DEP/


## distance between drug targets and INs (plasma specific network)
for reg in 'all' 'up' 'down'
do
    echo $reg
    mkdir -p $intermediate/REMAP/drug_proximity_plasma/IN/$reg

    python $src/REMAP/drug_proximity_forINs.py \
        -g $intermediate/intermediate_nodes/MPXV/plasma_proteome_from_Wang2022/intermediate_nodes_"$reg".txt \
        -d $result/REMAP/REMAP_pred30256.txt \
        -i $result/plasma_specific_network.txt \
        -o $intermediate/REMAP/drug_proximity_plasma/IN/$reg
done

## statistics of drug candidates for 3 protein sets
Rscript --vanilla $src/REMAP/drug_for_3sets.R

## network-based drug and side effect proximity (plasma specific network)
mkdir -p $intermediate/REMAP/drug_sideEffect_proximity_plasma

/share/home/dxj_tangkang/software/Python-3.9.12/python $src/REMAP/drug_sideEffect_proximity.py \
    -p $data/SIDER/side_effect_proteins.txt \
    -d $result/REMAP/REMAP_pred30256.txt \
    -l $result/REMAP/drug_proximity_plasma/drug_candidate.txt \
    -i $result/plasma_specific_network.txt \
    -o $intermediate/REMAP/drug_sideEffect_proximity_plasma/

## the final candidate drug list
Rscript --vanilla $src/REMAP/final_candidate_drugs.R























