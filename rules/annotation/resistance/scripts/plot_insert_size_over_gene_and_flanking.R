library(svglite)

data <- read.csv(snakemake@input[["insert_sizes"]], sep="\t", header=FALSE, stringsAsFactor=FALSE, col.names=c("Position", "InsertSize"))
bed <- read.csv(snakemake@input[["bed"]], sep="\t", header=FALSE, stringsAsFactor=FALSE, col.names=c("Reference", "Start", "End", "Strand"))

if(bed[, "Strand"] == -1) {title <- paste(snakemake@wildcards[["gene"]], "(negative strand)")
} else title <- paste(snakemake@wildcards[["gene"]], "(positive strand)")

increment <- seq(min(data[,1]), max(data[,1]), by=100)

sliding_standard_deviation <- c()
for (i in 1:(length(increment)-1)){
    sliding_standard_deviation <- c(sliding_standard_deviation, sd(data[data[,"Position"]>increment[i] & data[,"Position"]<increment[i+1], "InsertSize"]))
}

fileConn<-file(snakemake@output[["mean_insert_sizes"]])
writeLines(paste(snakemake@wildcards[["sample"]], as.integer(sd(data[,"InsertSize"]))), fileConn)
close(fileConn)

svglite(snakemake@output[["insert_sizes_plot"]])
plot(data[,"Position"], log(data[,"InsertSize"]), type="l", ylab="Standard deviation of reads insert sizes (log scale)", xlab=paste("Position on", snakemake@wildcards[["ref"]], "assembly"), main = title)
rect(bed[,"Start"] + snakemake@params[["up_down"]], 0, bed[, "End"] - snakemake@params[["up_down"]], log(max(data[,"InsertSize"])), col=rgb(96/256, 168/256, 98/256, alpha=0.5), border=NA)
graphics.off()
