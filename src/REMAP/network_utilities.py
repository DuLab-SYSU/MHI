import networkx
import random

def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree(): #.items(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = sorted(degree_to_nodes.keys())
    #values.sort()
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        #print i, low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins


def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes

def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False, seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                #nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20): # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [ random.choice(nodes) ]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values


# Imposes restrictions on func, func(a,b) should be always called as func(a, b, None) or func(a, b, dump_file=None)
def dumper(func):
    def wrapper(*args, **kwargs):
        import os
        def modify_kwargs():
	    # remove dump_file argument
            from copy import deepcopy
            kwargs_mod = deepcopy(kwargs) 
            del kwargs_mod["dump_file"]
            return kwargs_mod
	#print args, kwargs
        if "dump_file" in kwargs:
            #raise ValueError("dump_file argument not provided")
            dump_file = kwargs["dump_file"]
	    #kwargs_mod = modify_kwargs()
            kwargs_mod = kwargs
            args_mod = args
        else: # assuming the last argument as dump_file
            dump_file = args[-1]
	    #args_mod = args[:-1]
            kwargs_mod = kwargs
            args_mod = args
	#print args_mod, kwargs_mod
        if dump_file is not None:
            if os.path.exists(dump_file):
               val = cPickle.load(open(dump_file)) 
	       #print "loading", dump_file
            else:
		# remove dump_file argument
		#print "dumping", dump_file, kwargs_mod
                val = func(*args_mod, **kwargs_mod)
                if isinstance(val, types.GeneratorType): # newer versions of networkx returns iterators
                    val_new = {source: target_dict for source, target_dict in val}
                    val = val_new
                cPickle.dump(val, open(dump_file, 'w'))
        else:
	    #print "running", kwargs_mod
            val = func(*args_mod, **kwargs_mod)
        return val
    return wrapper

@dumper
def test(a, b, dump_file=None):
    print(a + b)
    return

@dumper
def get_shortest_paths(G, dump_file):
    return networkx.shortest_path(G)
    #if dump_file is not None:
    #	if os.path.exists(dump_file):
    #	   sp = cPickle.load(open(dump_file)) 
    #	else:
    #	    sp = networkx.shortest_path(G)
    #	    cPickle.dump(sp, open(dump_file, 'w'))
    #else:
    #	sp = networkx.shortest_path(G)
    #return sp

@dumper
def get_shortest_path_lengths(G, dump_file):
    return networkx.shortest_path_length(G)

@dumper
def get_shortest_path_lengths_subset(G, nodes_subset, dump_file):
    d = networkx.shortest_path_length(G)
    if nodes_subset is None:
        return d
    d_new = dict((node, {}) for node in nodes_subset)
    for node in nodes_subset:
        if node not in d:
            continue
        d_new[node] = dict((node2, d[node][node2]) for node2 in nodes_subset)
    return d_new

def get_shortest_path_between(G, source_id, target_id):
    return networkx.shortest_path(G, source_id, target_id)

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)

def get_all_shortest_paths_between(G, source_id, target_id):
    return networkx.all_shortest_paths(G, source_id, target_id)

def get_all_paths_from(G, source_id): 
    """
        get all paths from source node to all possible nodes in a dictionary
    """
    return networkx.single_source_dijkstra_path(G, source_id)

def get_all_paths_between(G, source_id, target_id, cutoff=None): 
    return networkx.all_simple_paths(G, source_id, target_id, cutoff)

def get_steiner_tree(G, terminals, sp=None):
    subgraph = NWSteiner.NWConnSteiner(G, terminals, shortestPath=sp)
    return subgraph

def get_heuristic_tree(G, terminals, distance=2, score_increment=0.2):
    # d: search radius, r: the expansion factor
    subgraph = Heuristic.listQuery(G, terminals, d=distance, r=score_increment) # d=1, r=0.2, scoreFun=Heuristic.sum_score
    return subgraph

def get_kwalk_tree(G, terminals):
    subgraph = kWalk.limkWalks(terminals, G) #, L=50, iteration=1)
    return subgraph

@dumper
def get_clustering_coefficient(g, dump_file):
    return networkx.clustering(g)


@dumper
def get_closeness_vitality(g, dump_file):
    return networkx.closeness_vitality(g)


def get_network_radius(g):
    return networkx.radius(g)


def get_network_degree_histogram(g):
    return networkx.degree_histogram(g)


@dumper
def get_node_betweenness(G, dump_file):
    return networkx.betweenness_centrality(G)

@dumper
def get_edge_betweenness(G, dump_file):
    return networkx.edge_betweenness_centrality(G)



def get_separation(network, sp, targets, seeds, distance, parameters={}, averaging_function=lambda x: numpy.mean(x)): 
    """
    Potential averaging functions
    val = numpy.min(x) 
    val = numpy.mean(x) 
    val = -numpy.log(numpy.mean([numpy.exp(-value-1) for value in x]))
    jorg / mmtom / avg distances: knn-x, shortest, kernel, mahalanobis, tom, mtom, ... 
    """
    if distance.startswith("jorg-"):
        distance = distance[len("jorg-"):] # closest was default before
        target_to_distance = get_source_to_average_target_distance(sp, seeds, seeds, distance = distance, parameters = parameters, target_mean_and_std = None, exclude_self=True)
        values = target_to_distance.values()
        d1 = averaging_function(values)
        target_to_distance = get_source_to_average_target_distance(sp, targets, targets, distance = distance, parameters = parameters, target_mean_and_std = None, exclude_self=True)
        values = target_to_distance.values()
        d2 = averaging_function(values)
        target_to_distance = get_source_to_average_target_distance(sp, targets, seeds, distance = distance, parameters = parameters, target_mean_and_std = None)
        values = target_to_distance.values()
        target_to_distance = get_source_to_average_target_distance(sp, seeds, targets, distance = distance, parameters = parameters, target_mean_and_std = None)
        values.extend(target_to_distance.values())
        d12 = averaging_function(values)
        val = d12 - (d1 + d2) / 2.0
    elif distance.startswith("mahalanobis-jorg-"):
        distance = distance[len("mahalanobis-jorg-"):] # mahalanobis is default before
        target_to_distance, d1 = get_source_to_average_target_distance(sp, targets, seeds, distance, parameters = parameters)
        values = target_to_distance.values()
        target_to_distance, d2 = get_source_to_average_target_distance(sp, seeds, targets, distance, parameters = parameters)
        values.extend(target_to_distance.values())
        d12 = averaging_function(values)
        val = d12 - (d1 + d2) / 2.0
    elif distance == "mahalanobis-pairwise":
        target_to_distance, center_d = get_source_to_average_target_distance(sp, targets, seeds, distance = "mahalanobis-shortest", parameters = parameters)
        values = target_to_distance.values()
        target_to_distance, center_d = get_source_to_average_target_distance(sp, seeds, targets, distance = "mahalanobis-shortest", parameters = parameters)
        values.extend(target_to_distance.values())
        val = averaging_function(values)
    elif distance == "center-pairwise":
        center_targets, center_values = get_center_of_subnetwork(sp, targets)
        center_seeds, center_values = get_center_of_subnetwork(sp, seeds)
        values = []
        for c_t in center_targets:
            for c_s in center_seeds:
                values.append(sp[c_t][c_s])
        val = averaging_function(values)
    elif distance == "closest-pairwise":
        values = []
        for geneid in targets:
            lengths = sp[geneid]
            inner_values = []
            for geneid_seed in seeds:
                if geneid == geneid_seed:
                    val = 0
                else:
                    val = lengths[geneid_seed]
                inner_values.append(val)
            values.append(min(inner_values))
        for geneid in seeds:
            lengths = sp[geneid]
            inner_values = []
            for geneid_target in targets:
                if geneid == geneid_target:
                    val = 0
                else:
                    val = lengths[geneid_target]
                inner_values.append(val)
            values.append(min(inner_values))
        val = averaging_function(values)
    elif distance == "shortest-pairwise":
        values = []
        for geneid in targets:
            lengths = sp[geneid]
            for geneid_seed in seeds:
                if geneid == geneid_seed:
                    val = 0
                else:
                    val = lengths[geneid_seed]
                values.append(val)
        val = averaging_function(values)
    elif distance == "dsd-pairwise" or distance == "communicability-pairwise":
        DSD, name_to_idx = sp
        values = []
        for geneid in targets:
            i = name_to_idx[geneid]
            for geneid_seed in seeds:
                if geneid == geneid_seed:
                    val = 0
                else:
                    j = name_to_idx[geneid_seed]
                    if distance.startswith("communicability"):
                        val = 1 - DSD[i, j]
                    else:
                        val = DSD[i, j]
                values.append(val)
        for geneid in seeds:
            i = name_to_idx[geneid]
            for geneid_target in targets:
                if geneid == geneid_target:
                    val = 0
                else:
                    j = name_to_idx[geneid_target]
                    val = DSD[i, j]
                values.append(val)
        val = averaging_function(values)
    elif distance == "kernel-pairwise":
        values = []
        for geneid in targets:
            lengths = sp[geneid]
            for geneid_seed in seeds:
                if geneid == geneid_seed:
                    val = 0
                else:
                    val = lengths[geneid_seed]
                values.append(val)
        val = -numpy.log(numpy.sum([numpy.exp(-value-1) for value in values])) / len(values)
    else:
        if distance.endswith("tom") or distance.endswith("tom-min"):
            target_to_distance = get_source_to_average_target_overlap(network, targets, seeds, distance)
        elif distance.startswith("mahalanobis"):
            target_to_distance, center_d = get_source_to_average_target_distance(sp, targets, seeds, distance, parameters = parameters)
        elif distance == "tsesolc":
            target_to_distance = get_source_to_average_target_distance(sp, seeds, targets, "closest", parameters = parameters)
        else:
            target_to_distance = get_source_to_average_target_distance(sp, targets, seeds, distance, parameters = parameters)
        values = target_to_distance.values()
        if distance.endswith("-min"):
            val = numpy.min(values)
        else:
            val = averaging_function(values)
    return val


