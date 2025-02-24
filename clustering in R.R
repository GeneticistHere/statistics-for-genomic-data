library(devtools)
library(Biobase)
library(dendextend)


con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()



edata = edata[rowMeans(edata) > 5000,]
edata = log2(edata + 1)

# By default calculates the distance between rows
dist1 = dist(t(edata))

## Look at distance matrix
colramp = colorRampPalette(c(6,"white",4))(9)
heatmap(as.matrix(dist1),col=colramp,Colv=NA,Rowv=NA)

hclust1 = hclust(dist1)
plot(hclust1)
plot(hclust1,hang=-1)


dend = as.dendrogram(hclust1)
dend = color_labels(hclust1,4,col=1:4)
plot(dend)


labels_colors(dend) = c(rep(1,10),rep(2,9))
plot(dend)


########kmean clustering
kmeans1 = kmeans(edata,centers=3)
names(kmeans1)
matplot(t(kmeans1$centers),col=1:3,type="l",lwd=3)
table(kmeans1$cluster)


heatmap(as.matrix(edata)[order(kmeans1$cluster),],col=colramp,Colv=NA,Rowv=NA)
kmeans2 = kmeans(edata,centers=3)
table(kmeans1$cluster,kmeans2$cluster)








