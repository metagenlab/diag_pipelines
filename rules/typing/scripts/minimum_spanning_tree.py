import pandas
import networkx 
import numpy
import matplotlib.pyplot as plt

matrix = pandas.read_csv(snakemake.input[0], sep="\t", index_col=0)

labels = {}
for i in range(len(matrix.index.values)):
    labels[i] = matrix.index.values[i]
    
a = numpy.array(matrix)
print(a)
k = networkx.from_numpy_matrix(a)

print(networkx.get_node_attributes(k, "label"))

T=networkx.minimum_spanning_tree(k)
edge_labels = networkx.get_edge_attributes(T,'weight')
pos=networkx.spring_layout(T)
plt.figure(figsize=(20,20))
networkx.draw_networkx_nodes(T, pos)
networkx.draw_networkx_edges(T, pos)
networkx.draw_networkx_edge_labels(T, pos, edge_labels)
networkx.draw_networkx_labels(T, pos, labels)
plt.savefig(snakemake.output[0])
