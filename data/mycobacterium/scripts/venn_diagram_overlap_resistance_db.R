library(VennDiagram)

rgi <- read.table("../db/rgi_annotated_full_2_0_0.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
rgi[,"MutatedAminoAcidOrNucleotide"] <- sapply(rgi[,"MutatedAminoAcidOrNucleotide"], toupper)
rgi[,"WildTypeAminoAcidOrNucleotide"] <- sapply(rgi[,"WildTypeAminoAcidOrNucleotide"], toupper)


mykrobe <- read.table("../db/mykrobe_annotated.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
mykrobe[,"MutatedAminoAcidOrNucleotide"] <- sapply(mykrobe[,"MutatedAminoAcidOrNucleotide"], toupper)
mykrobe[,"WildTypeAminoAcidOrNucleotide"] <- sapply(mykrobe[,"WildTypeAminoAcidOrNucleotide"], toupper)


miotto=read.table("../db/miotto_high_moderate_minimum_confidence_annotated.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
miotto[,"MutatedAminoAcidOrNucleotide"] <- sapply(miotto[,"MutatedAminoAcidOrNucleotide"], toupper)
miotto[,"WildTypeAminoAcidOrNucleotide"] <- sapply(miotto[,"WildTypeAminoAcidOrNucleotide"], toupper)


walker=read.table("../db/walker_resistant_annotated.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
walker[,"MutatedAminoAcidOrNucleotide"] <- sapply(walker[,"MutatedAminoAcidOrNucleotide"], toupper)
walker[,"WildTypeAminoAcidOrNucleotide"] <- sapply(walker[,"WildTypeAminoAcidOrNucleotide"], toupper)


genes <-  c("eis", "rpoB", "rpsL","gyrA", "gyrB", "katG", "pncA", "rrs", "tlyA", "gidB", "inhA", "embA", "embB", "embC", "embR", "ethA", "folC", "iniA", "iniB", "iniC", "kasA", "ndh", "ribD", "thyA", "ahpC", "mabA")

mutation_type <- "SNP"
for (gene in genes){
    rgi_res <- unique(sort(rgi[rgi[,"MutationType"]==mutation_type & rgi[,"Gene"]==gene,"PositionMTB"]))
    walker_res <-  unique(sort(walker[walker[,"MutationType"]==mutation_type & walker[,"Gene"]==gene,"PositionMTB"]))
    miotto_res <-  unique(sort(miotto[miotto[,"MutationType"]==mutation_type & miotto[,"Gene"]==gene,"PositionMTB"]))
    mykrobe_res <-  unique(sort(mykrobe[mykrobe[,"MutationType"]==mutation_type & mykrobe[,"Gene"]==gene,"PositionMTB"]))
    venn.diagram(list(CARD=rgi_res, Walker=walker_res, Miotto=miotto_res, Mykrobe=mykrobe_res), filename=paste("../venns/", gene, ".svg", sep=""), imagetype="svg", main=gene, fill=c("#cc545e","#64a860","#9970c1","#b98d3e"), cat.col=c("#cc545e","#64a860","#9970c1","#b98d3e"), width=10, height=10)
}
