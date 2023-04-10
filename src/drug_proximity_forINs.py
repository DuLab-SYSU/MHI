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
    description = 'Network-based separation of virus and pathway.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g', '--IN_gene', type=str, required=True, help='Intermediate nodes')
    parser.add_argument('-d', '--drug_target_file', type=str, required=True, help='Drug target filename')
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-o', '--out_dir', type=str, required=True, help='Statistic of proximity')
    #parser.add_argument('-o2', '--expected_lcc_size', type=str, required=True, help='Random expectation of LCC size')
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
    nodes_random = network_utilities.pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed) 
    return nodes_random

def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                val = network_utilities.get_shortest_path_length_between(network, node_from, node_to)
                values.append(val)
            d = min(values)
	    #print d,
            values_outer.append(d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append(val)
            d = min(values)
            values_outer.append(d)
    d = np.mean(values_outer)
    return d

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456, lengths=None, distance="closest"):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network 
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None # At least one of the node group not in network
    if distance != "closest":
        lengths = network_utilities.get_shortest_path_lengths(network, "temp_n%d_e%d.sif.pcl" % (len(nodes_network), network.number_of_edges()))
        d = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
    else:
        d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = network_utilities.get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = np.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        if distance != "closest":
            values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        else:
            values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    #pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    m, s = np.mean(values), np.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s) #(z, pval)

def get_drug_target(drug_target_file, network):
    drug_file = open(drug_target_file)
    drug_target = {}
    for row in itertools.islice(drug_file, 1, None):
        cols = row.strip().split('\t')
        drug = cols[2]
        if cols[5] in list(network.nodes()):
            drug_target.setdefault(drug, []).append(cols[5])
    return drug_target

def drug_proximity(network, drug_target, drug, nodes_to, out_dir):
    nodes_from = set(drug_target[drug])
    d, z, (m, s) = calculate_proximity(network, nodes_from, nodes_to)
    open(out_dir + '/' + drug + ".txt", 'a').write("%s\t%s\t%s\t%s\t%s\n" % ('drug', 'distance', 'z-score', 'distance_mean', 'distance_sd'))
    open(out_dir + '/' + drug + ".txt", 'a').write("%s\t%f\t%f\t%f\t%f\n" % (drug, d, z, m, s))
    return
    #return dict(zip([drug], [(d, z, m, s)]))

def run(args):
    G = create_network_from_first_two_columns(args.edge_list_file, delim = '\t')
    
    ins = []
    with open(args.IN_gene) as f:
        for row in f:
            ins.append(row.strip())

    ins = [i for i in ins if i in list(G.nodes())] # target genes

    start = time.time()
    drug_target = get_drug_target(args.drug_target_file, G)
    combs = list(itertools.product([G], [drug_target], [i for i in drug_target], [ins], [args.out_dir]))
    with mp.Pool(mp.cpu_count()-2) as pool:
        res = pool.starmap(drug_proximity, combs)
        pool.close()
        pool.join()
    end = time.time()
    print("The running time of this process is: {}h".format((end-start)/3600))

    #open(args.drug_distance, 'a').write("%s\t%s\t%s\t%s\t%s\n" % ('drug', 'distance', 'z-score', 'distance_mean', 'distance_sd'))
    #for i in res:
    #    for drug, stat in i.items():
    #        open(args.drug_distance, 'a').write("%s\t%f\t%f\t%f\t%f\n" % (drug, stat[0], stat[1], stat[2], stat[3]))
    

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))

