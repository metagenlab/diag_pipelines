#!/usr/bin/env python

def gc_coverage_plot(contigs_file,
                     contig_depth_table=False,
                     samtool_depth_file=False, 
                     blast_file=False,
                     column1=1,
                     column2=2,
                     main=False,
                     highlight=False,
                     taxonomy_file=False,
                     output_prefix=False):

    import os
    import shell_command
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from pandas import DataFrame
    import pandas
    import rpy2
    from rpy2.robjects import r
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    if not main:
        main = os.path.basename(contigs_file)
 
    out, err, code = shell_command.shell_command("infoseq -auto -only -Name -length -pgc %s > /tmp/gc.tab" % contigs_file)

    #print (out)
    #print (err)
    #print (code)

    if contig_depth_table:

        contig_depth = pandas.read_csv(contig_depth_table, sep='\t', names=["contig","depth"])
        #contig_depth = DataFrame(contig_depth, columns=['contig', 'depth'])
        #print (type(contig_depth["contig"]))
        #print (type(contig_depth))
        robjects.r.assign('contigs_depth', pandas2ri.py2ri(contig_depth))

    if taxonomy_file:
        with open(taxonomy_file, 'r') as f:
            contigs2taxon2count = {}

            for row in f:
                data = row.rstrip().split()
                contig = data[0]
                taxon = data[1]
                if contig not in contigs2taxon2count:
                    contigs2taxon2count[contig] = {}
                    contigs2taxon2count[contig][taxon] = 1
                else:
                    if taxon in contigs2taxon2count[contig]:
                        contigs2taxon2count[contig][taxon] +=1
                    else:
                        contigs2taxon2count[contig][taxon] = 1
        contig2label = []
        for contig in contigs2taxon2count:
            if len(contigs2taxon2count[contig]) > 1:
                # more than one taxon
                label = ''
                for taxon in contigs2taxon2count[contig]:
                    label+='%s (%s) /' % (taxon,
                                          contigs2taxon2count[contig][taxon])
                label = label[0:-2]
            else:
                label = list(contigs2taxon2count[contig].keys())[0]
            contig2label.append([contig,
                                 label])
        label2freq = {}
        for contig in contig2label:
            if contig[1] not in label2freq:
                label2freq[contig[1]] = 1
            else:
                label2freq[contig[1]] +=1
        for contig in contig2label:
            if label2freq[contig[1]] <=2:
                contig[1] = 'rare_taxon'

        df = DataFrame(contig2label, columns=['contig', 'label'])
        print (type(df["contig"]))
        print (type(df))
        #m = m.astype(float)
        robjects.r.assign('contig_labels', pandas2ri.py2ri(df))
    else:
        robjects.r.assign('contig_labels', False)

    if highlight:
        highlight_code = """
        gc_coverage_table$color <- rep(rgb(1, 0, 0,0.5), length(gc_coverage_table[,1]))
        highlight_table <- read.table("%s", header=FALSE)
        m <- match(highlight_table[,1], gc_coverage_table$Name)
        gc_coverage_subset <- gc_coverage_table[m,]
        print("subset")
        print(m)
        gc_coverage_table[m,]$color<-rgb(0, 0, 1,0.5)

        """ % highlight

        highlight_code2 = """

        m <- match(highlight_table[,1], gc_coverage_table_2m$Name)
        #print("subset m2")
        #print(m)
        gc_coverage_subset2 <- gc_coverage_table_2m[m,]

        """

    else:
        highlight_code = ''
        highlight_code2 = ''

    if not blast_file:
        robjects.r("""

        #library(Cairo)
        library(R.utils)
        library(ggplot2)



        
        if (exists("contigs_depth")==FALSE){
        
            if (isGzipped("%s")){
                #print('Gzipped file')
                all_depth <- read.table(gzfile('%s'), header=FALSE)
            }else{
                #print('Not Gzipped')
                all_depth <- read.table('%s', header=FALSE)
            }
      
            contigs_depth<- aggregate(all_depth["V3"],b=all_depth[c("V1")],FUN=median)
            colnames(contigs_depth) <- c('contig', 'depth')
        }
        #print(contigs_depth)
        #print(contig_labels)
        contigs_gc <- read.table("/tmp/gc.tab", header=TRUE)

        gc_coverage_table <-cbind(contigs_gc,coverage=contigs_depth[match(contigs_gc$Name, contigs_depth$contig),2])
        #w<-which(gc_coverage_table$Length >=1000)
        #gc_coverage_table <- gc_coverage_table[w,]

         cov_biggest <- gc_coverage_table[which(gc_coverage_table$Length==max(gc_coverage_table$Length)),4]
         #print('cov biggest:')
         #print(cov_biggest)
         w <- which(gc_coverage_table[,4]< (4*cov_biggest))
         gc_coverage_table_2m <- gc_coverage_table[w,]
    
        if (contig_labels != FALSE) {
            library(RColorBrewer)
            color_palette <- c('red', 'blue','green', brewer.pal(12,"Paired"), brewer.pal(12,"Set3"))
            m <- match(contig_labels$contig, gc_coverage_table$Name)

            gc_coverage_table$color <- rep("Unclassified", length(gc_coverage_table[,1]))
            gc_coverage_table$contig_alpha <- rep(0.5, length(gc_coverage_table[,1]))
        
            gc_coverage_table$color[m] <- as.character(contig_labels$label)
            gc_coverage_table$contig_alpha[m] <- rep(1,length(contig_labels$label))

            w<-which(gc_coverage_table$Length >=1000)
            gc_coverage_table <- gc_coverage_table[w,]
            #w2 <- which(gc_coverage_table$color != "Chlamydiae")
            #gc_coverage_table$contig_alpha[w2] <- 0.7

            svg("%sgc_cov_buble_test.svg", width = 12, height = 12)
            p6 <- ggplot(gc_coverage_table, aes(x = X.GC, y = coverage, size = Length, fill = color, colour = color, alpha = contig_alpha)) +
                    geom_point(shape = 21) +
                    ggtitle("Scaffold GC vs Depth") +
                    labs(x = "GC (%%)", y = "Sequencing depth") +
                    scale_size(range = c(1, 10))
            p6 <- p6 + scale_fill_manual(values=color_palette[0:length(unique(gc_coverage_table$color))])+ guides(color = guide_legend(override.aes = list(size=5)))
            p6 <- p6 + scale_colour_manual(values=color_palette[0:length(unique(gc_coverage_table$color))])
            #print (max(gc_coverage_table$Length))
            p6 <- p6 + scale_alpha_continuous(range=c(0.1, 1), limits=c(0.1,1)) #+ scale_alpha_continuous(range=c(0, max(gc_coverage_table$Length)), limits=c(0,max(gc_coverage_table$Length)))
            
            print(p6 + theme_bw())
            dev.off()

            gc_coverage_table_2m$color <- rep("Unclassified", length(gc_coverage_table_2m[,1]))
            gc_coverage_table_2m$contig_alpha <- rep(0.5, length(gc_coverage_table_2m[,1]))
        
            gc_coverage_table_2m$color[m] <- as.character(contig_labels$label)
            gc_coverage_table_2m$contig_alpha[m] <- rep(1,length(contig_labels$label))

            svg("%sgc_cov_buble_test_2m.svg", width = 12, height = 12)
            p6 <- ggplot(gc_coverage_table_2m, aes(x = X.GC, y = coverage, size = Length, fill = color, colour = color, alpha = contig_alpha)) +
                    geom_point(shape = 21) +
                    ggtitle("Scaffold GC vs Depth") +
                    labs(x = "GC (%%)", y = "Sequencing depth") +
                    scale_size(range = c(1, 10))
            p6 <- p6 + scale_fill_manual(values=color_palette[0:length(unique(gc_coverage_table$color))])+ guides(color = guide_legend(override.aes = list(size=5)))
            p6 <- p6 + scale_colour_manual(values=color_palette[0:length(unique(gc_coverage_table$color))])
            #print (max(gc_coverage_table$Length))
            p6 <- p6 + scale_alpha_continuous(range=c(0.1, 1), limits=c(0.1,1)) #+ scale_alpha_continuous(range=c(0, max(gc_coverage_table$Length)), limits=c(0,max(gc_coverage_table$Length)))
            
            print(p6 + theme_bw())
            dev.off()


            
        }else{
            #print('NO contig_labels')

        }
        
        write.table(gc_coverage_table, 'gc_coverage_table.tab', sep="\t", row.names=F)

    %s

     svg("%sgc_cov_buble.svg", width = 12, height = 12)
         symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3, ann=T,
                 bg=rgb(1, 0, 0,0.5), fg=rgb(1, 0, 0,0.5), main="%s", xlab="GC(%%)", ylab="Sequencing depth")
         if (any("gc_coverage_subset" %%in%% ls())) {
             symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3,
                     ann=T, bg=gc_coverage_table$color, fg=gc_coverage_table$color, add = TRUE)
             l <- gsub('(^[^_]+_[^_]+)_(.*)$', '\\\\1', gc_coverage_subset$Name)
             text(x=gc_coverage_subset[,3], y=gc_coverage_subset[,4], labels = l)
         }else{
            print ('a')
         }

         dev.off()




         %s

         svg("%sgc_cov_buble_2m.svg", width = 12, height = 12)
            symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2],
                    inches=1/3, ann=T, bg=rgb(1, 0, 0,0.5), fg=rgb(1, 0, 0,0.5), main="%s", xlab="GC(%%)", ylab="Sequencing depth")

            if (any("gc_coverage_subset" %%in%% ls())) {

                symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2],
                        inches=1/3, ann=T, bg=gc_coverage_table_2m$color, fg=gc_coverage_table_2m$color, add = TRUE)
                l <- gsub('(^[^_]+_[^_]+)_(.*)$', '\\\\1', gc_coverage_subset2$Name)
                text(x=gc_coverage_subset2[,3], y=gc_coverage_subset2[,4], labels = l)
            }else{
                print ('a')
            }

     dev.off()



                   """ % (samtool_depth_file,
                          samtool_depth_file,
                          samtool_depth_file,
                          output_prefix,
                          output_prefix,
                          highlight_code,
                          output_prefix,
                          main,
                          highlight_code2,
                          output_prefix,
                          main))
    else:

        robjects.r("""

        #library(Cairo)
        library(R.utils)

        if (isGzipped("%s")){
            #print('Gzipped file')
            all_depth <- read.table(gzfile('%s'), header=FALSE)
        }else{
            #print('Not Gzipped')
            all_depth <- read.table('%s', header=FALSE)
        }

        blast_file <- read.table("%s", header=FALSE, sep="\t")[,c(2,6)]
        contigs_depth<- aggregate(all_depth["V3"],b=all_depth[c("V1")],FUN=median)
        contigs_gc <- read.table("/tmp/gc.tab", header=TRUE)

        gc_coverage_table <-cbind(contigs_gc,coverage=contigs_depth[match(contigs_gc$Name, contigs_depth$V1),2])
        #w<-which(gc_coverage_table$Length >=1000)
        #gc_coverage_table <- gc_coverage_table[w,]

        gc_coverage_table$taxon <- blast_file[,2][match(gc_coverage_table$Name, blast_file[,1])]
        #print (is.na(gc_coverage_table$taxon))
        gc_coverage_table$taxon <- as.character(gc_coverage_table$taxon)
        gc_coverage_table$taxon[is.na(gc_coverage_table$taxon)] <- 'undefined'
        gc_coverage_table$taxon <- as.factor(gc_coverage_table$taxon)
        
        write.table(gc_coverage_table, 'gc_coverage_table.tab', sep="\t", row.names=F)

         svg("gc_cov_buble.svg", width = 12, height = 12,)
            symbols(x=gc_coverage_table[,3], y= gc_coverage_table[,4], circles=gc_coverage_table[,2], inches=1/3,
                    ann=F, bg=gc_coverage_table$taxon, fg=gc_coverage_table$taxon, main="%s", xlab="GC(%%)", ylab="Sequencing depth")
         dev.off()

         cov_biggest <- gc_coverage_table[which(gc_coverage_table$Length==max(gc_coverage_table$Length)),4]
         #print('cov biggest:')
         #print(cov_biggest)
         w <- which(gc_coverage_table[,4]< (4*cov_biggest))
         gc_coverage_table_2m <- gc_coverage_table[w,]

         svg("gc_cov_buble_2m.svg", width = 12, height = 12,)
            symbols(x=gc_coverage_table_2m[,3], y= gc_coverage_table_2m[,4], circles=gc_coverage_table_2m[,2],
                    inches=1/3, ann=F, bg=gc_coverage_table$taxon, fg=gc_coverage_table$taxon, main="%s", xlab="GC(%%)", ylab="Sequencing depth")
         dev.off()



                   """ % (samtool_depth_file,
                          samtool_depth_file,
                          samtool_depth_file,
                          blast_file,
                          main,
                          main))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_samtools_depth', type=str, help="input samtools depth files")
    parser.add_argument("-m", '--input_contigs', type=str, help="input contigs")
    parser.add_argument("-d", '--contig_depth_table', type=str, help="tab file with contig depth")
    parser.add_argument("-b", '--blast_file', type=str, help="file with blast columns")
    parser.add_argument("-1", '--column1', type=str, help="contig names column", default=2)
    parser.add_argument("-2", '--column2', type=str, help="classification column", default=6)
    parser.add_argument("-l", '--highlight', type=str, help="highlight some contigs listed in input file", default=False)
    parser.add_argument("-f", '--filter_size', type=int, help="filter contigs smaller than X", default=1000)
    parser.add_argument("-t", '--taxonomy', type=str, help="color contigs based on taxonomy", default=False)
    parser.add_argument("-o", '--output_prefix', type=str, help="output_prefix", default="./")
    
    args = parser.parse_args()

    gc_coverage_plot(args.input_contigs,
                     args.contig_depth_table,
                     args.input_samtools_depth,
                     args.blast_file,
                     args.column1,
                     args.column2,
                     False,
                     args.highlight,
                     taxonomy_file=args.taxonomy,
                     output_prefix=args.output_prefix)
