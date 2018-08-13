


# Read in libraries and set seed

rm(list=ls())
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR") #already installed on Aaron's macbookpro
# Load the GenVisR package
library("GenVisR");library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2)
set.seed(426)


setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP")




### Read in and organize main datasets FAP samples: clinical and mutation

        #Read in clinical information for adding additional information to waterfall plots
        sampleData <- read.table("https://docs.google.com/spreadsheets/d/e/2PACX-1vRd0tND0K2rcfq_C6pte4_W42vkuwsemkv_A_jnrm6mm5N72nxxw5Bezumy898XmtzEQlcxfrYW9j2E/pub?output=csv",sep = ",",header = T)
        
        ### Read in and reorganize the Somatic Mutation table. Prepare for Waterfall plot creation.
        file<-read.table(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified", header = T, sep = "\t",quote = "")
        #Select only statistically interesting columns. Rename the columns so that are compatible with sample information
        file_2<-select(file, chr, pos, id, ref, alt, geneName, geneLoc, functionalClass, somatic, germline,rep,sc,ends_with(match = "maf"), ends_with(match = "d"),-CADD_phred, -Polyphen2_HVAR_pred, -mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad)%>%
                filter(!(rep==1 & sc ==1))%>% #remove low grade mutations
                filter(is.na(germline))%>% #remove any germline mutations
                filter(!(rep==1 & geneLoc!="exonic")) #only allow mutations in repeat regions in exons.


        file_varreads<-file_2[,c(1:42)];rm(file_2) #selects variant descriptor columns and maf columns
        colnames(file_varreads)<-str_replace_all(colnames(file_varreads),pattern = "maf",replacement = "")
        
        
## filter out mutations that have both rep=1 and sc=1

### For loops to analyze different Exonic regions and Patient data.
patients <- c("JP","EP")
regions <- c("exonic","non-exonic")
for (pt in patients) {
        for (region in regions) {
                #choose only necessary column data. Replace "-" with "."
                #Chose patient EP to start with
                # clinical<-select(sampleData, "Sample_Name", "Stage", "Location_of_Sampled_Specimen", "Size_.mm.", "Size_SML.S..5.5.M..10.L.10.", "Grade", "PatientSample_ez",  "mean_coverageMetrics")%>%
                #         filter(grepl(paste(pt),PatientSample_ez))%>%
                #         mutate(sample=str_replace_all(PatientSample_ez, pattern = "-",replacement = "."))%>%
                #         select(-Sample_Name, Grade = Grade, Size = "Size_SML.S..5.5.M..10.L.10." , Location = Location_of_Sampled_Specimen , Stage = Stage)
                #Create the long melted version of the clinical table
                # clinical.melt<-melt(clinical,id.vars = "sample")
                # 
                # target <- c( "Size", "Stage")
                # clinical.melt<-filter(clinical.melt, variable %in% target)
                
                #select only variants per polyp which are known to be somatic
                file_varreads_melt<-gather(file_varreads,"SampleName", "MAF",13:42)%>%
                        filter(grepl(paste(pt),SampleName))%>%
                        filter(str_detect(somatic, paste0(SampleName,"\\","[")))%>% #filters for only somatic mutations per sample
                        filter(MAF>=0.05)
                
                # #Create new column for variant class - if in geneLoc exon, add the functional class.
                if (region=="exonic") {
                        #Create table with only exon functional variations (remove the non-exonic variants)
                        file_varreads_melt$variant_class <- file_varreads_melt$functionalClass
                        file_varreads_melt_2<-filter(file_varreads_melt, variant_class!="NA") 
                } else if (region=="non-exonic") {
                        #Create table without exon functional variations (remove the exonic variants)
                        file_varreads_melt$variant_class <- file_varreads_melt$geneLoc
                        file_varreads_melt_2<-filter(file_varreads_melt,variant_class!="NA")
                }
                
                
                #Making custom MAF plots per patient. make each varialbe a factor
                customMAF <- select(file_varreads_melt_2, "sample"="SampleName","gene"="geneName", variant_class )%>% #sample,gene,variant_class
                        mutate_all(funs(factor))
                
                #
                if (pt=="EP") {
                        #EP samples to plot - in order
                        samples_to_plot <- c("EP_Dec_NL", "EP.11",  "EP.31", "EP.57",  "EP.84","EP.26","EP.6", "EP.74B", "EP.88B", "EP.AdenoCa")        
                } else if (pt=="JP") {
                        #JP samples to plot - in order
                        samples_to_plot <- c("JP_Dec_NL.1", "JP69", "JP9",  "JP38","JP3", "JP31", "JP63", "JP34B", "JP61B", "JP6B", "JPAdenoCa")        
                }
                
                ########### For Making Waterfall Plots                
                # 
                #for exonic regions only
                main_layer <- theme_grey()
                clinical_layer <- theme_classic()+
                        theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10))

                # #make pdfs of waterfall plots. Exonic and Non-exonic need different mutations
                # pdf(file=paste0(pt,"_",region,"_variants_waterfall.pdf"), height=20, width=15)
                # if (region=="exonic") {
                #         waterfall(customMAF, fileType = "Custom",sampOrder = samples_to_plot, variant_class_order = c("nonsynonymous SNV", "stopgain", "synonymous SNV","unknown"), mainDropMut = T,mainXlabel = T,mainRecurCutoff = .03,plotMutBurden = T,mainPalette = c("purple","red", "green", "grey" ), mainLayer = main_layer, clinData = clinical.melt, clinLayer = clinical_layer, clinLegCol = 2)
                # } else if (region=="non-exonic") {
                #         customMAF_noexon<-filter(customMAF,!(variant_class=="exonic"))
                #         nonexonVariants <- c("intergenic", "intronic", "splicing", "upstream", "downstream", "UTR3", "UTR5", "ncRNA_exonic", "ncRNA_intronic")
                #         nonexonColors <- c("antiquewhite2", "brown", "darkolivegreen3", "deepskyblue2", "dodgerblue4", "darkgoldenrod3", "darkgoldenrod4", "gray32", "gray72")
                #         pdf(file=paste0(pt,"_",region,"_variants_waterfall.pdf"), height=20, width=15)
                #         waterfall(customMAF_noexon, fileType = "Custom", sampOrder = samples_to_plot,variant_class_order = nonexonVariants,mainDropMut = T,mainXlabel = T,mainRecurCutoff = .5,plotMutBurden = T,mainPalette = nonexonColors, mainLayer = main_layer, clinData = clinical.melt, clinLayer = clinical_layer, clinLegCol = 2)
                #         dev.off()
                # }
                
                #determine which genes are mutated more than once
                length(customMAF$gene) #Total number of mutations found in patient samples
                length(unique(customMAF$gene)) #number of distinct mutated genes in patient samples
                repeat_gene<-group_by(file_varreads_melt_2,geneName)%>%
                        arrange(geneName)%>%
                        count(geneName)%>% #counts number of times a gene is mutated (both identical and non-identical)
                        arrange(geneName)%>%
                        filter(n>1)%>%
                        arrange(geneName)
                        
                #creates table with shared mutated genes.
                table <- inner_join(x = file_varreads_melt_2, y = repeat_gene,"geneName")%>%
                        group_by(geneName)%>%
                        arrange(pos)%>%
                        arrange(-n)%>%
                        rename(times.gene.mutated=n)
                write.csv(table,file = paste0(pt,"_",region,"_sharemutatedgenes_ver2.csv"),append = F,col.names = T,sep=",",row.names = F)
                
                #creates table with shared point mutations
                table_2 <- add_count(table,pos, id, sort = T)%>%
                        rename(times.gene.mutated.identical=n)%>%
                        filter(times.gene.mutated.identical>1)%>%
                        arrange(-times.gene.mutated)%>%
                        arrange(-times.gene.mutated.identical)
                write.csv(table_2,file = paste0(pt,"_",region,"_sharemutatedgenes_identical_ver2.csv"),append = F,col.names = T,sep=",",row.names = F)
                
        }
}qqqq
        

## Pool all identically shared mutations together between EP and JP.
files=($(ls *sharemutatedgenes_identical_ver2*))
rm identicalMutations_ver2.csv #so that we dont continue to add new lines to the same file
cat  ${files[0]} ${files[1]} ${files[2]} ${files[3]} >> identicalMutations_ver2.csv

x<-read.csv(file = "identicalMutations_ver2.csv",sep = ",",header = T)%>%
        filter(!(chr=="chr"))%>%
        mutate_at(vars(times.gene.mutated),funs(as.numeric))%>%
        mutate_at(vars(times.gene.mutated.identical),funs(as.numeric))%>%
        mutate_at(vars(pos),funs(as.character))%>%
        mutate_at(vars(pos),funs(as.numeric))
        

bases=100
y<-select(x,-germline)%>%distinct()%>%arrange(-times.gene.mutated.identical)%>%arrange(-times.gene.mutated)#%>%
#         mutate(upstream = pos - bases)%>%
#         mutate(downstream = pos + bases)%>%
#         select(chr,upstream,downstream,pos,id,ref,alt,geneName,geneLoc,functionalClass,rep,sc,somatic,n,nn)#%>%
#         mutate(igv_prep = paste0("samtools view -b bam_file ","chr",chr,":",upstream,"-",downstream ))



write.csv(x = y,file = "identicalMutations_unique_ver2.csv",append = F,col.names = T,sep=",",row.names = F)
# less identicalMutations_ver2.csv | cut -f 1,2 -d , | grep -v "chr" | sort | uniq | wc -l #number of identical mutations

#Ruping annotated the file in the last column for repeat regions. I filtered them out

tempfile<-read.table(file = "identicalMutations_unique.tsv.rep", header = T, sep = "\t",quote = "") 
identical.nonrepeats<-tempfile[is.na(tempfile$test.rep.table.hg38.repeats_UCSC.withheader.gff.attribute),]
write.csv(x = identical.nonrepeats,file = "identicalMutations_unique_nonrepeats.csv",append = F,col.names = T,sep=",",row.names = F)
# maf <- select(file, contains('maf')) #select only "maf" minor allele frequency columns
# d <- select(file, ends_with('d'),-id,) #select only "d" depth columns
# num_var_reads <- maf*d #multiplying the two determines the number of reads which contained the variant.
# colnames(num_var_reads)<-str_replace_all(colnames(num_var_reads),pattern = "maf",replacement = "") #relabel the columns for the num_var_reads
# 
# file<-cbind(file,num_var_reads) #adds num_var_reads columns to file table
# file_varreads<-file[,c(1:10,71:100)] #selects variant descriptor columns and num_var_reads


# ######### Germline mutations #############
# 
# APC<-filter(file, chr==5,grepl(pattern = "APC",x = geneName))
# col_APC<-colnames(APC)
#         select(APC,col_APC,ends_with("maf"))
# 
# file_varreads_melt<-gather(file_varreads,"SampleName", "MAF",11:40)%>%
#         filter(grepl("EP",SampleName))%>%
#         filter(str_detect(somatic, paste0(SampleName,"\\","[")))%>% #filters for only somatic mutations per sample
#         filter(MAF>=0.05)