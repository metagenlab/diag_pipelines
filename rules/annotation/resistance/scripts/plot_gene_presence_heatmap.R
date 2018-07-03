library(ggplot2)
library(gridExtra)

rgi_files <- snakemake@input[["rgi_files"]]
species <- snakemake@params[["species"]]
mechanism <- snakemake@input[["rgi_mechanism_files"]]
print(mechanism)
sample_table <- data.frame(file=rgi_files, species=species,mechanism=mechanism, stringsAsFactors=FALSE)

for (i in 1:nrow(sample_table)){
  rgi_results_file <- sample_table[i,1]
  mechanism_file <-sample_table[i,3]
  species <- sample_table[i,2]

  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.csv(rgi_results_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    dataset$sample <- rep(unlist(strsplit(rgi_results_file, "/"))[2], length(dataset[,1]))
    dataset$species <- rep(species, length(dataset[,1]))
    resistance_mechanism <-  read.csv(mechanism_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    m <- match(dataset$Best_Hit_ARO, resistance_mechanism$Gene)
    dataset$mechanism <- resistance_mechanism$Resistance.Mechanism[m]


  }
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.csv(rgi_results_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    temp_dataset$sample <- rep(unlist(strsplit(rgi_results_file, "/"))[2], length(temp_dataset[,1]))
    temp_dataset$species <- rep(species, length(temp_dataset[,1]))
    resistance_mechanism <-  read.csv(mechanism_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    m <- match(temp_dataset$Best_Hit_ARO, resistance_mechanism$Gene)
    temp_dataset$mechanism <- resistance_mechanism$Resistance.Mechanism[m]

    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

# sort factors
u <- unique(dataset$Best_Hit_ARO)
u_sort <- u[rev(order(u))]
dataset$Best_Hit_ARO <- factor(dataset$Best_Hit_ARO, levels=u_sort)

# overview plot
pdf(snakemake@output[["rgi_plot"]], height=25,width=0.5*length(rgi_files))
p <- ggplot(data = dataset, aes(x = sample, y = Best_Hit_ARO)) + geom_raster(aes(fill = CUT_OFF))
p <- p + theme_grey(base_size = 10)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
print (p + coord_fixed(ratio=1))
dev.off()

###
# produce one plot per species
nr_species <- unique(sample_table$species)

for (species in nr_species){
    sub_dataset <- dataset[dataset$species==species,]
    sub_dataset$Best_Hit_ARO <- as.character(sub_dataset$Best_Hit_ARO)
    u <- unique(sub_dataset$Best_Hit_ARO)
    u_sort <- u[rev(order(u))]
    sub_dataset$Best_Hit_ARO <- factor(sub_dataset$Best_Hit_ARO, levels=u_sort)

    w <- length(unique(sub_dataset$species))*6
    h <- length(unique(sub_dataset$Best_Hit_ARO))/4
    pdf(paste('resistance/', species, '.pdf', sep=''), width=w, height=h)
        p <- ggplot(data = sub_dataset, aes(x = sample, y = Best_Hit_ARO)) + geom_tile(aes(fill = CUT_OFF), height = 0.9, width=0.9)
        p <- p + theme_grey(base_size = 10)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        print(p + coord_fixed(ratio=1)) #+ theme(strip.text.x = element_text(size=8, angle=75))
    dev.off()
}
