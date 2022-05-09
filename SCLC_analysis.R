rm (list= ls())
library(GEOquery)
library(limma)


################################################################
#Retrive gene expression data from GEO
gset <- getGEO("GSE149507", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL23270", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Sample details
pData(gset) #sample information
fData(gset) # sample annotation
exprs(gset)[1,]# expression data

# check whether data normalized or not
pData(gset)$data_processing[1]#data normalised

#  check the expression values
boxplot(exprs(gset),main="Distribution",outline=FALSE, las=2)
summary(exprs(gset))

################################################################
#Cleaning the samples
gset = na.exclude(gset)

#check the genes expression data whether contain NA values or not
gExprs <- exprs(gset)
apply(gExprs, 2, function(x) any(is.na(x)))# return FALSE if no NA values

################################################################
#Extract selected gene for project
library(readr)
library(dplyr)
features = fData(gset)
features = features %>% select(c("ID","ENTREZ_GENE_ID","Description"))  

#Selected gene 1 (ASCL1)
ASCL1_info = filter(features, ID =="429_at")
ASCL1_info 

#Selected gene 2 (NEUROD1)
NEUROD1_info = filter(features, ID =="4760_at")
NEUROD1_info

#Extract the expression values of selected genes
ASCL1 = gExprs["429_at",]
NEUROD1 = gExprs["4760_at",]
################################################################
# General expression data analysis
#Compute the mean and median for gene expression values
mean(ASCL1)
median(ASCL1) 

mean(NEUROD1)
median(NEUROD1)

#Checking for outliers using boxplot
boxplot(ASCL1, main="Boxplot of ASCL1 Expression Values",xlab ="ASCL1",
        ylab="Expression Values",ylim=c(0,15),outline=TRUE, las =1,col= "red")


boxplot(NEUROD1, main="Boxplot of NEUROD1 Expression Values",xlab ="NEUROD1",
        ylab="Expression Values",ylim=c(0,15),outline=TRUE, las =1,col= "blue")

#Checking the normality of data
qqnorm(ASCL1,main="QQ plot for ASCL1",pch =16,las = 1, col= "red")
qqline(ASCL1, col="darkred")
shapiro.test(ASCL1)

qqnorm(NEUROD1,main="QQ plot for NEUROD1",pch =16,las = 1, col= "blue")
qqline(ASCL1, col="blue4")
shapiro.test(NEUROD1)

################################################################
#Load the FASTA file for gene 1 and gene 2

#Download the fasta file from ncbi website 
#ASCL1 (https://www.ncbi.nlm.nih.gov/gene/?term=ASCL1)
#NEUROD1 (https://www.ncbi.nlm.nih.gov/gene/?term=NEUROD1)
#Rename the files to 'ASCL1' and 'NEUROD1'

library(seqinr)
ASCL1fas <- read.fasta(file="ASCL1.fasta")
NEUROD1fas <- read.fasta(file="NEUROD1.fasta")

#Assigning both genes into vectors
ASCL1seq <- ASCL1fas[[1]]
NEUROD1seq <- NEUROD1fas[[1]]

################################################################
#Pie chart of the nucleotide base composition for ASCL1 gene
summary(ASCL1seq)
slices <- c(36691, 26805, 26610, 38420)
lbls <- c("A", "C", "G", "T")
pct <- round (slices/sum(slices)*100)
lbls <- paste(lbls, pct) #to include percentages in the labels
lbls <- paste(lbls, "%", sep="") # adding % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)), 
    main="Nucleotide Base Composition of ASCL1")

#Pie chart of the nucleotide base composition for NEUROD1 gene
summary(NEUROD1seq)
slices <- c(5529, 4103, 4121, 5470)
lbls <- c("A", "C", "G", "T")
pct <- round (slices/sum(slices)*100)
lbls <- paste(lbls, pct) #to include percentages in the labels
lbls <- paste(lbls, "%", sep="") # adding % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)), 
    main="Nucleotide Base Composition of NEUROD1")

legend ("bottomright", c("A", "C", "G", "T"), 
        cex = 0.5, 
        fill = c("red", "green", "light blue", "purple"), 
        title = "Type of nucleotide bases")

################################################################
#GC Content
GC(ASCL1seq)
GC(NEUROD1seq)

ASCL1gc <- (14+2)*100/(57)   #GCpercentage
NEUROD1gc <- (17+6)*100/(65) #GCpercentage
d <-c(ASCL1gc,NEUROD1gc)
matricgc <- matrix(d, byrow=TRUE, nrow = 1, ncol = 2)
colnames(matricgc) <- c("ASCL1","NEUROD1")
rownames(matricgc) <- c("GC Content (%)")
matricgc

barplot(main= "GC Content Percentage for ASCL1 and NEUROD1",
        c(ASCL1gc,NEUROD1gc), names.arg = c("ASCL1","NEUROD1"),
        xlab = "Gene", ylab = "GC Content (%)", las = 1,
        col = c("#cc9be7","#60e5c0"))
