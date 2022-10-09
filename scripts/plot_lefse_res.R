#!/usr/bin/env Rscript

library(ggplot2)

options(bitmapType='cairo')


read_data <- function(data){

    data <- read.table(data, sep="\t", header=TRUE)
    colnames(data) <- c("Taxonomy", "Log", "Group", "LDA", "Pvalue")
 
    #修改物种名称
    n <- 0
    taxs <- c()    
    for (i in data[,1]){
        n <- n+1
        temp <- strsplit(i, "__")
        np <- length(temp[[1]])
        taxs <- c(taxs, temp[[1]][np])

    }
    data[,1] <- taxs 
    return(data)

}


plot_lefse_res <- function(data, prefix, top=20, alphorder=FALSE){

    data <- read_data(data)
    if(alphorder){
        data <- data[order(data$Group, data$Taxonomy, data$LDA), ] #按照分组、名字和LDA值进行排序
    }else {
        data <- data[order(data$Group, data$LDA), ]
    }

    groups <- unique(data$Group)
    gnum <- length(groups)

    data[data$Group==groups[gnum],4] = 0 - data[data$Group==groups[gnum], 4] #将最后一个组转化为负数
    
    k <- FALSE
    for(i in groups){
        if(k){
            temp <- rbind(temp, head(data[data$Group==i,], top))
        }else{
            temp <- head(data[data$Group==i,], top)
            k <- TRUE
        }
    }
    data <- temp
    num <- length(data[,1])
    if(num <= 2){
        height <- 5
    }else {
        height <- 5 + (num-2)*0.3
    }
    if(height <=50){
       height <= 50
    }
    data <- data[!duplicated(data$Taxonomy),] #去除重复
    data$Taxonomy <- factor(data$Taxonomy, levels=as.character(data$Taxonomy))
 
    p <- ggplot(data, aes(x=Taxonomy, y=LDA, fill=Group)) +
        geom_bar(stat='identity', colour='black', width=0.9, position=position_dodge(0.7))+
        xlab("") + ylab("LDA") + coord_flip() +
        theme_bw() +
        theme(axis.text.y=element_blank(), axis.ticks=element_blank(),
          axis.text.x=element_text(face='bold'), axis.title.x=element_text(face='bold')) +
        theme(panel.border=element_blank()) + 
        geom_text(aes(y=ifelse(data$LDA >0,-0.1,0.1), label=Taxonomy), fontface=4, size=3, hjust=ifelse(data$LDA>0,1,0))
    
    ggsave(paste(prefix, "lefse.png", sep="."), units="cm", width=14, height=height)
    ggsave(paste(prefix, "lefse.pdf", sep="."),  units="cm", width=14, height=height)
}


add_help_args <- function(args){

    if(length(args) <=1) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Plot lefse.\n")
        cat("Example1:plot_lefse_res.R lefse.res prefix\n")
        cat("Example2:plot_lefse_res.R lefse.res prefix 30\n")
        cat("The number of species on display is 30.\n")
        cat("Example4:plot_lefse_res.R lefse.res prefix 30 TRUE\n")
        cat("Sort by species letter.\n")
        cat("Input file format:\n")
        cat("\tTaxonomy\tLogarithm value\tGroup\tLDA-value\tP-value\n")
        cat("Escherichia_coli\t4.4132\tSM\t4.0341\t0.02836\n")
        quit()
    }
}

args <- commandArgs(trailingOnly=TRUE)
add_help_args(args)

if(length(args) ==2) {
    plot_lefse_res(args[1], args[2])
}else if(length(args) ==3){
    plot_lefse_res(args[1], args[2], args[3])
}else{
    plot_lefse_res(args[1], args[2], args[3], args[4])
}
