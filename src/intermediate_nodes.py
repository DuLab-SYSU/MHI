#!/usr/bin/python

import networkx as nx
import pandas as pd
import itertools
import sys
from tqdm import tqdm
import multiprocessing as mp
from multiprocessing import Pool
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests
import argparse

def get_parser():
    description = 'Identification of intermediate genes through shortest paths'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--virus_host_ppi', type=str, required=True, help='Virus-host PPI')
    parser.add_argument('-d', '--deg_file', type=str, required=True, help='Differential expressed genes')
    parser.add_argument('-r', '--regulation', type=str, required=True, help='Up and down regulation')
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-o', '--intermediated_nodes', type=str, required=True, help='Intermediated nodes')
    parser.add_argument('-o2', '--shortest_path', type=str, required=True, help='Interset shortest paths')
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
    
def my_func(G, source):
    res = nx.single_source_dijkstra_path(G, source)
    res = {k:v for k, v in res.items() if len(v)>2}
    res_dict = {source: res}
    return res_dict


def starmap_parallel(func, G, process_num=10):
    nodes = list(G)
    tasks = list(zip([G]*len(nodes), nodes))
    with Pool(process_num) as pool:
        results = pool.starmap(func, tqdm(tasks, total=len(tasks)))

    final = {}
    for v in results:
        final.update(v)
    return final



def get_shortest_path(all_paths_dict, sources, targets):
    print("Get all shortest paths")
    print('...')
    all_paths = []
    for key_source, path_dict in all_paths_dict.items():
        for key_target, path in path_dict.items():
            if path:
                all_paths.append(path)
    print('Finish all shortest path')
    print('')

    intrest_paths = []
    print("Get intrest paths")
    print('...')
    for source_i, target_i in itertools.product(sources, targets):
        path = all_paths_dict.get(source_i, {}).get(target_i)
        if not path:
            continue
        intrest_paths.append(path)
    print("Finsh intrest shortest path")
    return all_paths, intrest_paths


def get_in_node(path_set):
    return set([node for path in path_set for node in path[1:-1]])
    
def get_deg_node(deg_file, reg = 'all'):
    nodes = []
    with open(deg_file) as f:
        for row in f:
            if row.startswith("#"):
                pass
            else:
                col = row.strip().split('\t')
                if reg != 'all':
                    if col[1] == reg:
                        nodes.append(col[0])
                else:
                    nodes.append(col[0])
    return nodes

def single_hypergeo_test(node_id, all_set, intrest_set, M, N):
    n = sum(list(map(lambda x: 1 if node_id in x else 0, all_set)))
    k = sum(list(map(lambda x: 1 if node_id in x else 0, intrest_set)))
    p_value = stats.hypergeom.pmf(k=k, M=M, n=n, N=N)
    return node_id, p_value


def hypergeo_test_parallel(func, nodes, all_set, intrest_set, M, N, process_num=10):
    num_ = len(nodes)
    tasks = list(zip(nodes, [all_set]*num_, [intrest_set]*num_, [M]*num_, [N]*num_))
    with Pool(process_num) as pool:
        results = pool.starmap(func, tqdm(tasks, total=len(tasks)))

    final = []
    for v in results:
        final.append(v)
    return final
    

def run(args):
    G = create_network_from_first_two_columns(args.edge_list_file, delim = '\t')

    vhi = pd.read_csv(args.virus_host_ppi, sep='\t')
    vtg = set(vhi['Human Protein'].tolist())
    vtg = [i for i in vtg if i in list(G.nodes())]  # source genes
    
    deg = get_deg_node(args.deg_file, reg = args.regulation)
    deg = [i for i in deg if i in list(G.nodes())]
    
    #print(nx.info(G))
    #print('')
    
    # get_shortest_path_parallel
    print('All pairs shortest paths')
    all_paths_dict = starmap_parallel(my_func, G, process_num=mp.cpu_count())
    print('')
    all_shortest_paths, intreset_shortest_paths = get_shortest_path(all_paths_dict, vtg, deg)
    
    
    M = len(all_shortest_paths)
    N = len(intreset_shortest_paths)
    print('')
    print("Number of all paths: ", M)
    print("Number of intrest paths: ", N)

    in_nodes = list(get_in_node(intreset_shortest_paths))
    print("Number of intrest nodes: ", len(in_nodes))

    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('Significant Test: \n')
    pvalues = []
    for node in tqdm(in_nodes, total=len(in_nodes)):
        p = single_hypergeo_test(node, all_shortest_paths, intreset_shortest_paths, M, N)
        pvalues.append(p)
 
    #pvalues = hypergeo_test_parallel(single_hypergeo_test, in_nodes, all_shortest_paths, intreset_shortest_paths, M, N, process_num=mp.cpu_count())
    print('')
    print('Internal Node\tP value')
    for k, v in pvalues:
        print(k, '\t', v)

    index_, pvalues_ = zip(*pvalues)

    condi = multipletests(pvalues_, 0.01, 'fdr_bh')[0]
    edge_sig = np.array(index_)[condi]
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('Multi-correcting\n\n')
    print('Significant Internal Nodes: ')
    print(' '.join(edge_sig.tolist()))
    
    with open(args.intermediated_nodes, 'w') as f:
        f.write('\n'.join(edge_sig.tolist()))

    with open(args.shortest_path, 'w') as g:
        for x in intreset_shortest_paths:
            g.write(','.join(x) + '\n')
    
if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
