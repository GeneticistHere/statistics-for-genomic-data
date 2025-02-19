library(devtools)
library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()
par(mfrow=c(1,1))
hist(rnorm(1000),col=2)
hist(edata[,1],col=2,breaks=100)
hist(log(edata[,1]),col=2,breaks=100)
min(log(edata))
min(log(edata[,1] + 1))
hist(log(edata[,1] + 1),breaks=100,col=2)
hist(log2(edata[,1] + 1),breaks=100,col=2)
hist(log2(edata[,1] + 1),breaks=100,col=2,xlim=c(1,15),ylim=c(0,400))
hist(rowSums(edata==0),col=2)
low_genes = rowMeans(edata) < 5
table(low_genes)
filt_edata = filter(edata,!low_genes)
dim(filt_edata)
low_genes2 = rowMedians(as.matrix(edata)) < 5
table(low_genes2,low_genes)
# Instead of filter(), use direct subsetting
filt_edata = edata[!low_genes, ]
filt_edata2 = edata[!low_genes2, ]
dim(filt_edata2)
hist(log2(filt_edata[,1] + 1),col=2)







