library(svglite)

data <- read.csv(snakemake@input[["coverage_txt"]], sep="\t", header=FALSE, stringsAsFactor=FALSE, col.names = c("Reference", "Start", "End", "Strand", "Position", "Coverage"))

strand <- data[1, "Strand"]

if(strand == -1) {title <- paste(snakemake@wildcards[["gene"]], "(negative strand)")
} else title <- paste(snakemake@wildcards[["gene"]], "(positive strand)")

svglite(snakemake@output[["coverage_plot"]])
plot(data[,"Start"]+data[,"Position"], data[,"Coverage"], type="l", xlab=paste("Position on", snakemake@wildcards[["ref"]], "assembly"), ylab="Coverage", main = title)
rect(data[1, "Start"]+snakemake@params[["up_down"]], 0, data[1, "End"]-snakemake@params[["up_down"]], max(data[,"Coverage"]), col=rgb(96/256,168/256,98/256, alpha=0.5), border=NA)
graphics.off()

