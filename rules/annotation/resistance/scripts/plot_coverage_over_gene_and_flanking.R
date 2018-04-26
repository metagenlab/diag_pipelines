library(svglite)

data <- read.csv(snakemake@input[["coverage_txt"]], sep="\t", header=FALSE, stringsAsFactor=FALSE)

strand <- data[1, 4]

print(head(data))

if(strand == -1) {title <- paste(snakemake@wildcards[["gene"]], "(negative strand)")
} else title <- paste(snakemake@wildcards[["gene"]], "(positive strand)")

svg(snakemake@output[["coverage_plot"]])
plot(data$V2+data$V5, data$V6, type="l", xlab=paste("Position on", snakemake@wildcards[["ref"]], "assembly"), ylab="Coverage", main = title)
rect(data$V2[1]+snakemake@params[["up_down"]], 0, data$V3[1]-snakemake@params[["up_down"]], max(data$V6), col=rgb(96/256,168/256,98/256, alpha=0.5), border=NA)
graphics.off()

