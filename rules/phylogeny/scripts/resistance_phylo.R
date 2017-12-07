library(adephylo)
library(phylobase)

tree <- read.tree(snakemake.input[1])

resistance_data <- read.table(snakemake.input[2], sep=",", header=T, rownames=1)

resistance_data[resistance_data=="R"] <- 1
resistance_data[resistance_data=="S"] <- 0

data <- data.matrix(resistance_data)

p <- phylo4d(tree, data)

pdf(file=snakemake.output[1])

table.phylo4d(p,symbol="colors",scale=F,col=c("white","black"),legend=F,grid=T,show.tip.label=T,ratio.tree=0.15,cex.symbol=1,pch=15,box=F, cex.label=.75, center=FALSE, col.tip=c(rep("red", 31), rep("black", 29)), show.node.label=F)

dev.off()
