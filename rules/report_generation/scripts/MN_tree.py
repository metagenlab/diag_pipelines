

import numpy
import networkx
import pandas
import json
import itertools


def find_clusters(G, node_list):
    comb = itertools.combinations(node_list, 2)
    groups = []
    for i in comb:
        if i[0] == i[1]:
            continue
        try:
            print (G[i[0]][i[1]])
        except:
            if len(groups)==0:
                no_match=True
            else:
                no_match = False
            for n, group in enumerate(groups):
                if i[0] in group and i[1] in group:
                    continue
                elif i[0] in group and i[1] not in group:
                    groups[n].append(i[1])
                elif i[1] in group and i[0] not in group:
                    groups[n].append(i[1])
                else:
                    no_match=True
            if no_match:
                groups.append([i[0], i[1]])
    return groups

def get_spaced_colors(n):

    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

def merge_group_nodes(G, groups, node_list):
    for group in groups:
        median_dico = {}
        for node in node_list:
            if node in group:
                continue
            data = []
            for member in group:
                data.append(G[member][node]['weight'])
            m = numpy.median(data)
            median_dico[node] = m
        mapping = {group[0]: '\n'.join(group)}
        G = networkx.relabel_nodes(G, mapping)
        for i in group[1:len(group)]:
            G.remove_node(i)
        for i in node_list:
            if i in group:
                continue
            G['\n'.join(group)][i]['weight'] = median_dico[i]
    return G


# this function is used to convert networkx to Cytoscape.js JSON format
# returns string of JSON
def convert2cytoscapeJSON(G, node2st=False):

    if node2st:
        mlst_list = list(set(node2st.values()))
        mlst2color = dict(zip(mlst_list, get_spaced_colors(len(mlst_list))))
        mlst2color['-'] = 'white'
    print(mlst2color)
    # load all nodes into nodes array
    final = {}
    final["nodes"] = []
    final["edges"] = []
    for node in G.nodes():
        nx = {}
        nx["data"] = {}
        nx["data"]["id"] = node
        nx["data"]["label"] = node
        final["nodes"].append(nx.copy())
        if node2st:
            try:
                st = node2st[int(node.split('\n')[0])]
            except:
                st = '-'
            nx["data"]["color"] = mlst2color[st]
            nx["data"]["MLST"] = ST
        else:
            nx["data"]["color"] = '#8A8A8A'
    #load all edges to edges array
    for edge in G.edges(data=True):
        #print(edge)
        nx = {}
        nx["data"]={}
        nx["data"]["id"]=edge[0]+edge[1]
        nx["data"]["strength"] = edge[2]["weight"]
        if edge[2]["weight"] < 50:
            nx["data"]["color"] = '#f92411'
            nx["data"]["width"] = 2
        else:
            nx["data"]["color"] = '#8A8A8A'
            nx["data"]["width"] = 1
        nx["data"]["source"]=edge[0]
        nx["data"]["target"]=edge[1]
        final["edges"].append(nx)
    return json.dumps(final)


def get_MN_tree(dist_matrix):

    m = pandas.read_csv(dist_matrix, delimiter='\t', header=0, index_col=0)

    nodes = list(m.index)

    G = networkx.from_pandas_adjacency(m)  # networkx.from_numpy_matrix(A)

    groups = find_clusters(G, nodes)

    G = merge_group_nodes(G, groups, nodes)

    T=networkx.minimum_spanning_tree(G)

    return T