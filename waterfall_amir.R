


# Read in libraries and set seed

rm(list=ls())
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR") #already installed on Aaron's macbookpro
# Load the GenVisR package
library("GenVisR");library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2)
set.seed(426)


setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/AnnotationHive_Amir/")


#Making custom MAF plots per patient. make each varialbe a factor
                x<- read.csv(file = "customMAF_amir.csv",header = T,sep = ",")
                customMAF <- select(x, sample,gene, variant_class )%>% #sample,gene,variant_class
                        mutate_all(funs(factor))
                nonsynVariants <- c("nonsynonymous SNV")
                #make the plot into a pdf
                pdf(file=paste0("EP_nonsynonymous_Germline_AnnoHive.pdf"), height=20, width=15)
                waterfall(customMAF, fileType = "Custom", mainDropMut = T,mainXlabel = T,mainRecurCutoff = .5,plotMutBurden = T, variant_class_order = nonsynVariants)
                dev.off()
                
                