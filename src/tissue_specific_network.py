#!/usr/bin/python

import pandas as pd
import numpy as np
import networkx as nx
import random
import network_utilities
import sys, argparse

def get_parser():
    description = 'Calculate the tissue-specific network'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-g', '--tissue_gene', type=str, required=True, help='Tissue expressed genes')
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-o', '--subnetwork', type=str, required=True, help='Specific network')
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

   
def get_nodes(nodes_file):
    nodes = []
    with open(nodes_file) as f:
        for row in f:
            nodes.append(row.strip())
    return nodes


def run(args):
    G = create_network_from_first_two_columns(args.edge_list_file, delim = '\t')
    nodes = get_nodes(args.tissue_gene)
    
    sub_networks = G.subgraph(nodes)
    result_list = [c for c in sorted(nx.connected_components(sub_networks), key=len, reverse=True)]
    max_S = nx.to_pandas_edgelist(G.subgraph(result_list[0])) # the largest connected component
    max_S.to_csv(path_or_buf=args.subnetwork,
             index = 0, sep = '\t', header=0)
    

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))

