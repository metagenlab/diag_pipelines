data <- list()
for (i in snakemake@input){
    ref = gsub("alignments/flat_distance_", "", gsub("\\.txt", "", i))
    data[[i]] <- read.table(i, sep="\t", header=FALSE, row.names=1)
    }

data <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = 0, all = TRUE), data)

log1 <- replace(data[,3], data[,3]<1, 0.5)
log2 <- replace(data[,2], data[,2]<1, 0.5)

pdf(snakemake@output[[1]])
plot(log1~log2, log = "xy")
dev.off()
