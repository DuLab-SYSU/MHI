#!/usr/bin/python

import pandas as pd
import numpy as np
import networkx as nx
import random
import network_utilities
import sys, argparse

def get_parser():
    description = 'Calculate the viral module.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--virus_host_ppi', type=str, required=True, help='Virus-host PPI')
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-o1', '--lcc_size', type=str, required=True, help='Statistic of LCC')
    parser.add_argument('-o2', '--expected_lcc_size', type=str, required=True, help='Random expectation of LCC size')
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

def get_connected_components(G, return_as_graph_list=True):
    result_list = []

    if return_as_graph_list:
        result_list = nx.connected_component_subgraphs(G)
    else:
        result_list = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]

    return result_list

def lcc_size(network, nodes):
    network_sub = network.subgraph(nodes)
    component_nodes = get_connected_components(network_sub, False)
    d = len(component_nodes[0])
    return d

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100, degree_aware=True, seed=None):
    if bins is None:
        bins = network_utilities.get_degree_binning(network, min_bin_size)
    nodes_random = network_utilities.pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed) 
    return nodes_random

def calculate_lcc_significance(network, nodes, nodes_random=None, bins=None, n_random=1000, min_bin_size=100, seed=452456):
    if bins is None and nodes_random is None:
        bins = network_utilities.get_degree_binning(network, min_bin_size)
    
    random.seed(seed)
    if nodes_random is None:
        network_nodes = list(network.nodes())
        nodes_random = get_random_nodes(nodes, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    
    #network_sub = network.subgraph(nodes)
    #component_nodes = get_connected_components(network_sub, False)
    #d = len(component_nodes[0])
    d = lcc_size(network, nodes)
    
    values = np.empty(len(nodes_random))
    for i, nodes in enumerate(nodes_random):
        #network_sub = network.subgraph(nodes)
        #component_nodes = get_connected_components(network_sub, False)[0]
        #values[i] = len(component_nodes)
        values[i] = lcc_size(network, nodes)
    m, s = np.mean(values), np.std(values)
    if s == 0:
        z = 0
    else:
        z = (d - m) / s
    return d, z, (m, s), values

def run(args):
    G = create_network_from_first_two_columns(args.edge_list_file, delim = '\t')

    vhi = pd.read_csv(args.virus_host_ppi, sep='\t')
    vtg = set(vhi['Human Protein'].tolist())
    vtg = [i for i in vtg if i in list(G.nodes())]

    d, z, (m, s), values = calculate_lcc_significance(G, vtg)

    #print(d)
    with open(args.lcc_size, 'a') as out_f1, open(args.expected_lcc_size, 'a') as out_f2:
        out_f1.write("%s\t%s\t%s\t%s\n" % ('size', 'z-score', 'size_mean', 'size_sd'))
        out_f1.write("%d\t%f\t%f\t%f\n" % (d, z, m, s))
        for i in values:
            out_f2.write("%d\n" % (i))
    

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))

