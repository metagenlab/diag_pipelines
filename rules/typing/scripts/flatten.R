res <- read.table(snakemake@input[[1]])
res <- res[!rownames(res) == snakemake@wildcards[["ref"]], !colnames(res) == snakemake@wildcards[["ref"]]]
res <- res[order(rownames(res)),order(colnames(res))]

tmp <- c()
for (i in 1:dim(res)[1]){
    tmp <- c(tmp, setNames(res[i,], paste(colnames(res)[i], names(res[i,]), sep="_")))
}

tmp <- unlist(tmp)
itself <- paste(colnames(res), colnames(res), sep="_")
tmp <- tmp[!names(tmp) %in% itself]

write.table(tmp, snakemake@output[[1]], quote=FALSE, sep="\t", col.names=FALSE)
