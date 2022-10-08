#!/usr/bin/env Rscript

library("ggplot2")
library("ggpubr")
library("ggsci")
#library("ggprism")
#library("argparse")

options(bitmapType='cairo') #关闭服务器与界面的互动响应


read_abundance <- function(data, species){

   data <- read.delim(data, row.names=1, sep="\t", head=TRUE, check.names=FALSE)

   for (i in rownames(data)){
       if (grepl(species, i)){
           data <- data[i, ]
           rownames(data) <- c(species)
           break
       }
   }
   return(data)
}


get_compare_group <- function(groups){

    groups <- unique(groups)

    compare <- list()
    n <- 0
    for (i in 1:(length(groups)-1)) {
        for (j in (i+1):length(groups)){
            n <- n+1
            compare[[n]] <- c(groups[i], groups[j])
        }
    }

    return(compare)

}


plot_abundance_test_bar <- function(data, group, species){
    
    data <- read_abundance(data, species)
    data <- t(data)
    sample_id <- row.names(data)
    group <- read.delim(group, sep="\t", stringsAsFactors=FALSE, header=TRUE)
    colnames(group) <- c("sample", "group")
    rownames(group) <- group$sample
    group <- group[sample_id, ]
    group$sample <- data[,1]

    compare <- get_compare_group(group$group)
    maxylim <- max(data[,1])+max(data[,1])*0.2
    grouplen <- length(unique(group$group))
    if(grouplen <= 2){
        width <- 4
    }else {
        width <- 4+grouplen-2
        maxylim <- max(data[,1])+max(data[,1])*0.2*grouplen
    }


    p <- ggplot(group, aes(x=group, y=sample, fill=group)) +
        xlab("") +ylab("Relative abundance (%)") +
        #stat_summary(geom="col", fun=mean) +
        geom_bar(stat="summary", fun=mean, position="dodge",width=0.8) +
        stat_summary(geom="errorbar", fun=mean,
            fun.min= function(x) mean(x) - sd(x),
            fun.max= function(x) mean(x) + sd(x),
            width=0.3) +
        theme(panel.background=element_blank(), axis.line=element_line(), legend.position="none") +
        stat_compare_means(comparisons=compare, method="t.test", label="p.signif") +
        scale_y_continuous(expand=c(0, 0)) +
        coord_cartesian(ylim=c(0, maxylim))
    #print(compare) 
    ggsave(paste(species, "abundance_test_bar.png", sep="."), units="cm", width=width, height=6)
    ggsave(paste(species, "abundance_test_bar.pdf", sep="."), units="cm", width=width, height=6)

    colnames(group) <- c("abundance", "group")
    write.table(data.frame("sample"=rownames(group),group), paste(species, ".abundance.tsv", sep=""), row.names=FALSE, sep="\t", quote=FALSE, na ="")
}

add_help_args <- function(args){

    if(length(args) != 3) {
        cat("Version: v1.0.0\n")
        cat("Author:Boya Xu, Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:anosim analysis.\n")
        cat("Example:plot_abundance_test_bar.R abundance_species.xls group.list species\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9\n")
        quit()
    }
}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

plot_abundance_test_bar(args[1], args[2], args[3])
