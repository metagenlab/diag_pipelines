

library(igraph, warn.conflicts=FALSE)

set.seed(1)

sample_sts <- read.csv(snakemake@input[[2]], sep="\t", header=FALSE, row.names=1, stringsAsFactors=FALSE)
refs_sts <- read.csv(snakemake@input[[3]], sep="\t", header=FALSE, row.names=1, stringsAsFactors=FALSE)
all_sts <- rbind(sample_sts, refs_sts)

reference_name <-  read.csv(snakemake@input[[4]], sep="\t", header=FALSE, stringsAsFactors=FALSE)[1,1]

ref_convinient_name <- paste(reference_name, " (", snakemake@wildcards[["ref"]], ")", sep="")

print(ref_convinient_name)

#If the ST types files are empty : no MLST schema for the species under consideration
if (all(all_sts[,1]=="-")){
    sts_exist <- FALSE
} else {
    sts_exist <- TRUE
}

matrix <- as.matrix(read.csv(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1))
#We set the zero distance (clones) to a small value, otherwise no links exist between clones
matrix[matrix==0] <- 0.01

nb_vertex <- dim(matrix)[1]

graph <- graph.adjacency(matrix, mode="undirected", weighted=TRUE, diag=FALSE)
graph <- set_vertex_attr(graph, "name", value=gsub("X", "", vertex_attr(graph, "name")))


if (sts_exist){
    sts <- all_sts[vertex_attr(graph, "name"), 2]
    sts[sts == "-"] <- NA
    graph <- set_vertex_attr(graph, "ST", value=sts)
}

graph <- set_vertex_attr(graph, "name", value=gsub(snakemake@wildcards[["ref"]], ref_convinient_name, vertex_attr(graph, "name")))

mapping <- 1:nb_vertex

#We create the groups of clones
s1 <- clusters(subgraph.edges(graph, E(graph)[E(graph)$weight==0.01], delete.vertices=FALSE))



to_be_deleted <- c()
mapping <- 1:nb_vertex

#We prepare the vectors for grouping the clones and deleting the useless vertices
for (i in 1:s1$no){
    if (s1$csize[i] > 1){
        mapping[which(s1$membership==i)]=which(s1$membership==i)[1]
        to_be_deleted <- c(to_be_deleted, tail(which(s1$membership==i), -1))
    }
}

nb_vertex <- nb_vertex - length(to_be_deleted)


#Merging clones
graph <- contract(graph, mapping, vertex.attr.comb = "concat")
#Deleting useless vertices
graph <- delete_vertices(graph, to_be_deleted)


#Creation of the minimum spanning tree graph
mst_graph <- mst(graph, algorithm="prim")



#We keep all distances that are less than a particular threshold
for (i in E(graph)[E(graph)$weight<snakemake@params["threshold"]]){
    vert <- get.edges(graph, i)
    if (!are_adjacent(mst_graph, vert[,1], vert[,2]) && (vert[,1]!=vert[,2])){
        mst_graph <- add_edges(mst_graph, vert, weight=edge_attr(graph, "weight")[i])
    }
}

mst_graph <- set_edge_attr(mst_graph, "label", value=edge_attr(mst_graph, "weight"))

#We transform the weight so that smaller distances weight more. I tried an inverse square root and inverse log transformations, inverse log appears better.
mst_graph <- set_edge_attr(mst_graph, "weight", value=1/log(0.5+edge_attr(mst_graph, "weight")))

vertices_sizes <- sapply(vertex_attr(mst_graph, "name"), length, simplify=TRUE)*5
#Vertices sizes are proportional to the number of sample they represent


#Concatanating the names of the vertices to represent every clonal samples
mst_graph <- set_vertex_attr(mst_graph, "name", value=lapply(vertex_attr(mst_graph, "name"), paste, collapse="\n"))
mst_graph <- set_vertex_attr(mst_graph, "ST", value=lapply(vertex_attr(mst_graph, "ST"), unique))

vertices_colors = rep(NA, nb_vertex)



if (sts_exist){
    uniq_sts = unique(sts[!is.na(sts)])
    number_cols <- length(uniq_sts)
    colors <- rainbow(number_cols)
    for (i in 1:length(uniq_sts)){
    vertices_colors[vertex_attr(mst_graph, "ST")==uniq_sts[i]] <- colors[i]
    }
}


svg(snakemake@output[[1]], height=10, width=10)
if (sts_exist) {
    par(mar=c(0,8,0,0))
}



plot(mst_graph, vertex.size=vertices_sizes, edge.label.dist=35, vertex.label.dist=0, node.label.cex=0.5, vertex.color=vertices_colors, vertex.label=vertex_attr(mst_graph, "label"), vertex.label.color=rep("black", nb_vertex), vertex.label.family="sans", edge.label.family="sans")

if (sts_exist){
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("left", legend=paste("ST ", uniq_sts), fill=colors, xpd=TRUE)
}

dev.off()
