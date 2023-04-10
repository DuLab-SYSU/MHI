#!/usr/bin/python

import pandas as pd
import numpy as np
import networkx as nx
import random
import itertools
import time
import multiprocessing as mp
import network_utilities
import sys, argparse

def get_parser():
    description = 'Proximity algorithm for drug repurposing.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--virus_host_ppi', type=str, required=True, help='Virus-host PPI')
    parser.add_argument('-p', '--pathway_gene_file', type=str, required=True, help='Pathway gene filename')
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-o', '--pathway_separation', type=str, required=True, help='Pathway separation')
    parser.add_argument('-o2', '--expected_separation', type=str, required=True, help='Statistic of proximity')
    return parser


def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = nx.DiGraph()
    else:
        g = nx.Graph()
    return g

def create_network_from_first_two_columns(network_file, delim = None):
    g = create_graph()
    for line in open(network_file):
        id1, id2 = line.strip().split(delim)[:2]
        g.add_edge(id1, id2)
    return g

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        bins = network_utilities.get_degree_binning(network, min_bin_size)
    nodes_random = network_utilities.pick_random_nodes_matching_selected(network, bins, set(nodes), n_random, degree_aware, seed=seed) 
    return nodes_random

def get_separation(network, nodes_from, nodes_to, lengths=None):
    dAA = np.mean(get_separation_within_set(network, nodes_from, lengths))
    dBB = np.mean(get_separation_within_set(network, nodes_to, lengths))
    dAB = np.mean(get_separation_between_sets(network, nodes_from, nodes_to, lengths))
    d = dAB - (dAA + dBB) / 2.0
    return d


def get_separation_between_sets(network, nodes_from, nodes_to, lengths=None):
    """
    Calculate dAB in separation metric proposed by Menche et al. 2015
    """
    values = []
    target_to_values = {}
    source_to_values = {}
    for source_id in nodes_from:
        for target_id in nodes_to:
            if lengths is not None:
                d = lengths[source_id][target_id] 
            else:
                d = network_utilities.get_shortest_path_length_between(network, source_id, target_id)
            source_to_values.setdefault(source_id, []).append(d)
            target_to_values.setdefault(target_id, []).append(d)
    # Distances to closest node in nodes_to (B) from nodes_from (A)
    for source_id in nodes_from:
        inner_values = source_to_values[source_id]
        values.append(np.min(inner_values))
    # Distances to closest node in nodes_from (A) from nodes_to (B)
    for target_id in nodes_to:
        inner_values = target_to_values[target_id]
        values.append(np.min(inner_values))
    return values


def get_separation_within_set(network, nodes_from, lengths=None):
    """
    Calculate dAA or dBB in separation metric proposed by Menche et al. 2015
    """
    if len(nodes_from) == 1:
        return [ 0 ]
    values = []
    # Distance to closest node within the set (A or B)
    for source_id in nodes_from:
        inner_values = []
        for target_id in nodes_from:
            if source_id == target_id:
                continue
            if lengths is not None:
                d = lengths[source_id][target_id] 
            else:
                d = network_utilities.get_shortest_path_length_between(network, source_id, target_id)
            inner_values.append(d)
        values.append(np.min(inner_values))
    return values

def calculate_observed_separation_proximity(network, nodes_from, nodes_to, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set(network.nodes())
    if len(set(nodes_from) & nodes_network) == 0 or len(set(nodes_to) & nodes_network) == 0:
        return None # At least one of the node group not in network
    d = get_separation(network, set(nodes_from), set(nodes_to), lengths)
    return d

def calculate_random_separation_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set(network.nodes())
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = network_utilities.get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    
    random_values_list = zip(nodes_from_random, nodes_to_random)
    combs = [(network,i,j) for i,j in random_values_list]
    with mp.Pool(mp.cpu_count()-1) as pool:
        res = pool.starmap(get_separation, combs)
        pool.close()
        pool.join()
    m, s = np.mean(res), np.std(res)
    return m, s

def calculate_z_score(d, m, s):
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return z

def get_pathway_gene(pathway_gene_file, network):
    pathway_file = open(pathway_gene_file)
    pathway_gene = {}
    for row in itertools.islice(pathway_file, 1, None):
        cols = row.strip().split('\t')
        pathway = cols[1]
        if cols[2] in list(network.nodes()):
            pathway_gene.setdefault(pathway, []).append(cols[2])
    return pathway_gene

def pathway_separation_proximity(network, pathway_gene, pathway, nodes_to):
    """
    nodes_to here means virally-targeted genes
    """
    nodes_from = pathway_gene[pathway]
    d = calculate_observed_separation_proximity(network, nodes_from, nodes_to)
    return dict(zip([pathway], [(d)]))

def run(args):
    G = create_network_from_first_two_columns(args.edge_list_file, delim = '\t')

    vhi = pd.read_csv(args.virus_host_ppi, sep='\t')
    vtg = set(vhi['Human Protein'].tolist())
    vtg = [i for i in vtg if i in list(G.nodes())]

    #start = time.time()
    pathway_gene = get_pathway_gene(args.pathway_gene_file, G)
    combs = list(itertools.product([G], [pathway_gene], [i for i in pathway_gene], [vtg]))
    with mp.Pool(mp.cpu_count()-1) as pool:
        res = pool.starmap(pathway_separation_proximity, combs)
        pool.close()
        pool.join()
    #end = time.time()
    #print("The running time of this process is: {}h".format((end-start)/3600))

    
    open(args.pathway_separation, 'a').write("%s\t%s\n" % ('pathway', 'separation'))
    open(args.expected_separation, 'a').write("%s\t%s\t%s\t%s\t%s\n" % ('pathway', 'separation', 'z-score', 'separation_mean', 'separation_sd'))
    for i in res:
        for pathway, separation in i.items():
            open(args.pathway_separation, 'a').write("%s\t%f\n" % (pathway, separation))
            if separation < 0:
                nodes_from = pathway_gene[pathway]
                m, s = calculate_random_separation_proximity(G, nodes_from, vtg)
                z = calculate_z_score(separation, m, s)
                open(args.expected_separation, 'a').write("%s\t%f\t%f\t%f\t%f\n" % (pathway, separation, z, m, s))


if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))


