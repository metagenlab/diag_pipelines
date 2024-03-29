library(ggplot2)
library(gridExtra)
library(grid)
library(egg)
library(svglite)

rgi_files <- snakemake@input[["rgi_files"]]
species <- snakemake@params[["species"]]
sample_table <- data.frame(file=rgi_files, species=species, stringsAsFactors=FALSE)

for (i in 1:nrow(sample_table)){
  rgi_results_file <- sample_table[i,1]
  mechanism_file <-sample_table[i,3]
  species <- sample_table[i,2]

  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.csv(rgi_results_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    dataset$sample <- rep(unlist(strsplit(rgi_results_file, "/"))[2], length(dataset[,1]))
    dataset$species <- rep(species, length(dataset[,1]))
  }
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.csv(rgi_results_file, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    temp_dataset$sample <- rep(unlist(strsplit(rgi_results_file, "/"))[2], length(temp_dataset[,1]))
    temp_dataset$species <- rep(species, length(temp_dataset[,1]))
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

# sort factors
u <- unique(dataset$Best_hit)
u_sort <- u[rev(order(u))]
dataset$Best_hit <- factor(dataset$Best_hit, levels=u_sort)

# overview plot
# prepare one plot/resistance mechanism
plot_list <- list()
for (i in 1:length(unique(dataset$Mechanism))){
    resistance_mechanism <- unique(dataset$Mechanism)[i]
    mechanism_subset <- dataset[dataset$Mechanism==resistance_mechanism,]
    p <- ggplot(data = mechanism_subset, aes(x = sample, y = Best_hit)) + geom_tile(aes(fill = Cut_Off), height = 0.9, width=0.9)
    p <- p + theme_grey(base_size = 10)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p <- p + coord_fixed(ratio=1) + theme(legend.position="none") + ggtitle(resistance_mechanism) #+ theme(strip.text.x = element_text(size=8, angle=75))
    p <- p + theme(axis.title.x=element_blank())
    p <- p + theme(axis.title.y=element_blank())
    plot_list[[i]] <- p
}
# with of the plot should take into accounts:
# 1) the number of trains
# 2) the length of labels
label_length <- lapply(as.character(dataset$Best_hit), nchar)
longest_label <- max(unlist(label_length))
svglite(snakemake@output[["rgi_plot"]], height=25,width=longest_label/15 + (0.5*length(rgi_files)))
#p <- ggplot(data = dataset, aes(x = sample, y = Best_hit)) + geom_raster(aes(fill = CUT_OFF))
#p <- p + theme_grey(base_size = 10)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#print (p + coord_fixed(ratio=1))
ggarrange(plots=plot_list, ncol = 1, newpage = FALSE)
dev.off()

###
# produce one plot per species
nr_species <- unique(sample_table$species)

for (species in nr_species){
    sub_dataset <- dataset[dataset$species==species,]
    # sort factors
    sub_dataset$Best_hit <- as.character(sub_dataset$Best_hit)
    u <- unique(sub_dataset$Best_hit)
    u_sort <- u[rev(order(u))]
    sub_dataset$Best_hit <- factor(sub_dataset$Best_hit, levels=u_sort)
    # prepare one plot/resistance mechanism
    plot_list <- list()
    for (i in 1:length(unique(sub_dataset$Mechanism))){
        print("Mechanism")
        if(i == 0)
            next
        resistance_mechanism <- unique(sub_dataset$Mechanism)[i]
        mechanism_subset <- sub_dataset[sub_dataset$Mechanism==resistance_mechanism,]
        p <- ggplot(data = mechanism_subset, aes(x = sample, y = Best_hit)) + geom_tile(aes(fill = Cut_Off), height = 0.9, width=0.9)
        p <- p + theme_grey(base_size = 10)  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        p <- p + coord_fixed(ratio=1) + theme(legend.position="none") + ggtitle(resistance_mechanism) #+ theme(strip.text.x = element_text(size=8, angle=75))
        p <- p +   theme(axis.title.x=element_blank())
        p <- p +   theme(axis.title.y=element_blank())
        plot_list[[i]] <- p
    }

    # plot multiplot
    w <- length(unique(sub_dataset$species))*9
    h <- length(unique(sub_dataset$Best_hit))/2
    svglite(paste('report/resistance/', species, '.svg', sep=''), width=w, height=h)
        ggarrange(plots=plot_list, ncol = 1, newpage = FALSE)
    dev.off()
}
