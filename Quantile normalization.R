library(devtools)
library(Biobase)
library(preprocessCore)
library(ggplot2)
con =urlggplot2con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

#Show distributions for log2 counts for several samples
#Here we show density plots for the first 20 samples
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 3, ]
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}

p1<-plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}

#Quantile normalization
norm_edata = normalize.quantiles(as.matrix(edata))
p2<-plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.20))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}


#Matching distributions leaves variability
p3<-plot(norm_edata[1,],col=as.numeric(pdata$study))
svd1 = svd(norm_edata - rowMeans(norm_edata))
p4<-plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata$study))
#################Saving plots
# Save Plot 1
png("density_plots_before_normalization.png", width = 1600, height = 1600, res = 250)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30), main="Density Plots (Before Normalization)")
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}
dev.off()

# Save Plot 2
png("density_plots_after_normalization.png", width = 1600, height = 1600, res = 250)
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.20), main="Density Plots (After Quantile Normalization)")
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}
dev.off()

# Save Plot 3
png("scatter_plot_normalized_data.png", width = 1600, height = 1600, res = 250)
plot(norm_edata[1,],col=as.numeric(pdata$study), main="Scatter Plot of Normalized Data")
dev.off()

# Save Plot 4
png("pca_plot.png", width = 1600, height = 1600, res = 250)
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata$study), main="PCA Plot")
dev.off()







