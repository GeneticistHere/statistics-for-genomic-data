library(gplots)
library(ggplot2)
library(devtools)
library(Biobase)
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(RSkittleBrewer)
# Make the colors pretty
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)



con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
#to save this data for later use
save(bm, file = "bodymap_eset.xlsx")

pdata=pData(bm)
edata=exprs(bm)
fdata = fData(bm)
ls()


head(pdata)
dim(pdata)
#Tables for factor/character variables
table(pdata$gender)
table(pdata$gender,pdata$race)


#looking for NA values in edata
summary(edata)
# Use option useNA to include NA's in table
table(pdata$age,useNA="ifany")
# is.na checks for NA values
table(is.na(pdata$age))
# Check for other common missing names
sum(pdata$age==" ")
# Check genomic data for NAs
sum(is.na(edata))
# Make the distribution of NA's by genes
gene_na = rowSums(is.na(edata))
table(gene_na)
# Make the distribution of NA's by samples
sample_na = rowSums(is.na(edata))
table(sample_na)


dim(fdata)
dim(pdata)
dim(edata)
boxplot(log2(edata+1),col=2,range=0)

par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)


plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3)


qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)

mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)


edata = as.data.frame(edata)
filt_edata = filter(edata,rowMeans(edata) > 1)
dim(filt_edata)
boxplot(as.matrix(log2(filt_edata+1)),col=2)

aeid = as.character(fdata[,1])
chr = AnnotationDbi::select(org.Hs.eg.db,keys=aeid,keytype="ENSEMBL",columns="CHR")
head(chr)
dim(chr)
dim(edata)
# Take non-duplicated chromsomes
chr = chr[!duplicated(chr[,1]),]
# Confirm that the annotation still is in the right order
all(chr[,1] == rownames(edata))
# Select the chromosome Y samples
edatay = dplyr::filter(edata,chr$CHR=="Y")

# Males have Y chromsome expression as expected
boxplot(colSums(edatay) ~ pdata$gender)
points(colSums(edatay) ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender),
       pch=19)


ematrix = as.matrix(edata)[rowMeans(edata) > 10000,]
heatmap(ematrix)


colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(ematrix,col=colramp)


heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)
heatmap.2(ematrix,col=colramp,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none")


