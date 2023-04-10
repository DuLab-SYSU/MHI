#!/usr/bin/python

import pandas as pd
import networkx as nx
import glob, os
import re
import sys, argparse

def get_parser():
    description = 'Calculate the node properties for network.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-o', '--node_info_file', type=str, required=True, help='Node info filename')
    return parser

def run(args):
    f_index = re.split("[_.]", args.edge_list_file)[-2]
    G = nx.read_weighted_edgelist(args.edge_list_file, delimiter='\t')
    nodes = sorted(list(G.nodes()))

    """For each node"""
    nodes = sorted(list(G.nodes()))

    Degree = nx.degree(G)  # Degree
    k = [Degree[node] for node in nodes]

    Betweenness_centrality = nx.betweenness_centrality(G)  # Betweenness
    BC = [Betweenness_centrality[node] for node in nodes]

    Eigenvector_centrality = nx.eigenvector_centrality(G)  # Eigenvector centrality
    x = [Eigenvector_centrality[node] for node in nodes]

    Clustering_coefficient = nx.clustering(G) # local clustering coefficient
    C = [Clustering_coefficient[node] for node in nodes]

    Assortativity = nx.average_neighbor_degree(G) # The average neighborhood degree of a node
    NC = [Assortativity[node] for node in nodes]

    Closeness_centrality = nx.closeness_centrality(G) # reciprocal of the average shortest path distance
    SP = [Closeness_centrality[node] for node in nodes]

    node_info = pd.DataFrame({'Gene': nodes,
                              'Degree': k,
                              'Betweenness_centrality': BC,
                              'Eigenvector_centrality': x,
                              'Clustering_coefficient': C,
                              'Assortativity': NC,
                              'Closeness_centrality': SP,
                              'File_index': [f_index] * len(nodes)})
    
    node_info.to_csv(path_or_buf=args.node_info_file,
                     index = 0, 
                     sep = '\t')
    
if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
