library(igraph)


matrix <- as.matrix(read.csv(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1))

graph <- mst(graph_from_adjacency_matrix(matrix, mode="undirected", weight=TRUE), algorithm="prim")
graph <- set_edge_attr(graph, "label", value=edge_attr(graph, "weight"))
graph <- set_edge_attr(graph, "weight", value=1/sqrt(edge_attr(graph, "weight")+0.01))

print(dim(matrix))
graph <- set_vertex_attr(graph, "label", value=gsub("X", "", colnames(matrix)))

pdf(snakemake@output[[1]], height=10, width=10)
plot(mst(graph, algorithm="prim"), elab=TRUE, vertex.size=rep(3, dim(matrix)[1]), edge.label.dist=35, vertex.label.dist=3, node.label.cex=0.5, vertex.color=rep(NULL, 11), label.color=rep("black", 11))
dev.off()
