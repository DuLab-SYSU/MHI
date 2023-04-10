import re
import pandas as pd
import networkx as nx
import os

# The relationship of child and parent PSI MI (2022-05-18)
### The child MI and its recent parent MI
mi_relation = {}
child = []
with open(r'PSI-MI/mi.owl') as f:
    for row in f:
        if row.startswith('id: MI'):
            col = row.strip().split(" ")
            child.append(col[1])
        elif row.startswith('is_a: MI'):
            term = re.search('MI:[0-9]+', row.strip()).group()
            mi_relation.setdefault(child[-1], []).append(term)
            
child = []
parent = []
for ch in mi_relation.keys():
    for pa in mi_relation[ch]:
        child.append(ch)
        parent.append(pa)
        
mi_relation = pd.DataFrame({'Child': child,
                            'Parent': parent})

# The grandparent and parent of a child MI
### Find all the tip nodes of a given source
def find_ancestor(G, child):
    ancestors = []
    parents = nx.predecessor(G, child, cutoff = 1)
    p = list(parents.keys())
    if len(p) == 1:
        ancestor = p[0]
        ancestors.append(p[0])
    else:
        for node in p[1:]:
            ancestors += find_ancestor(G, node)
    return list(set(ancestors))

G = nx.from_pandas_edgelist(mi_relation, source='Child', target='Parent', 
                            edge_attr=None, create_using=nx.DiGraph())
### 'molecular interaction', 'interaction detection method', and 'experimental interaction detection' were not included
G.remove_nodes_from(['MI:0000','MI:0001','MI:0045'])

### The MI retained
mi_retain = []
for mi in list(G.nodes()):
    ancestors = find_ancestor(G, mi)
    for a in ancestors:
        if a in ['MI:0401', 'MI:0013', 'MI:0254', 'MI:0428', 'MI:1088', 'MI:0255', 'MI:0090']:
            mi_retain.append(mi)

mi_retain = set(mi_retain)
mi_remove = [i for i in list(G.nodes()) if i not in mi_retain]
G.remove_nodes_from(mi_remove)

G2 = G.copy()
G2.remove_nodes_from(['MI:0401', 'MI:0013', 'MI:0254', 
                     'MI:0428', 'MI:1088', 'MI:0255', 'MI:0090'])

### The functions
def find_grandparent(child):
    if child in list(G.nodes()):
        return find_ancestor(G, child)
    else:
        return "-"

def find_parent(child):
    if child in list(G2.nodes()):
        return find_ancestor(G2, child)
    else:
        return "-"

def miLevel(mi, a):
    return nx.shortest_path_length(G, mi, a, weight = None) + 1

