library(svglite)

data <- read.csv(snakemake@input[["insert_sizes"]], sep="\t", header=FALSE, stringsAsFactor=FALSE)
bed <- read.csv(snakemake@input[["bed"]], sep="\t", header=FALSE, stringsAsFactor=FALSE)

print(bed)
strand <- bed$V4

if(strand == -1) {title <- paste(snakemake@wildcards[["gene"]], "(negative strand)")
} else title <- paste(snakemake@wildcards[["gene"]], "(positive strand)")

increment <- seq(min(data[,1]), max(data[,1]), by=100)

sliding <- c()
for (i in 1:(length(increment)-1)){
    sliding <- c(sliding, sd(data[data$V1>increment[i] & data$V1<increment[i+1],]$V2))
}



fileConn<-file(snakemake@output[["mean_insert_sizes"]])
writeLines(paste(snakemake@wildcards[["sample"]], as.integer(mean(sliding))), fileConn)
close(fileConn)

svg(snakemake@output[["insert_sizes_plot"]])
plot(head(increment, -1), sliding, type="l", ylab="Standard deviation of reads insert sizes", xlab=paste("Position on", snakemake@wildcards[["ref"]], "assembly"), main = title)
rect(bed$V2+snakemake@params[["up_down"]], 0, bed$V3[1]-snakemake@params[["up_down"]], max(data$V2), col=rgb(96/256,168/256,98/256, alpha=0.5), border=NA)
graphics.off()


