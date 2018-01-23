

l <- read.table(snakemake@input[[1]])

min_value = Inf
index = 1
for (i in 1:500){
    s = min(which((cumsum(l[-(1:i),2])/sum(l[-1,2])>0.95)==TRUE))
    if (s <= min_value){
        min_value = s
        index = i
    }
}

write(index, snakemake@output[[1]])
write(min_value+index, snakemake@output[[2]])

