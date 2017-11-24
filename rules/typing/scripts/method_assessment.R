data <- read

print(snakemake@input)

tmp <- c()
for (i in 1:dim(res)[1]){
    tmp <- c(tmp, setNames(res[i,], paste(colnames(res)[i], names(res[i,]), sep="_")))
    }
        
    
