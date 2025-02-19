# Required Libraries
library(devtools)
library(Biobase)

# Install required packages
install.packages(c("devtools"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase"))

# Load gene expression data from URL
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)

# Create data objects
mp = montpick.eset
pdata = pData(mp)  # phenotype data
edata = as.data.frame(exprs(mp))  # expression data
fdata = fData(mp)  # feature data

# Filter for genes with mean expression > 100
edata = edata[rowMeans(edata) > 100, ]

# Log transform data
edata = log2(edata + 1)

# Center data by subtracting row means
edata_centered = edata - rowMeans(edata)

# Calculate SVD
svd1 = svd(edata_centered)

# Plot singular values
png("singular_values.png", width = 1600, height = 1600, res = 250)
plot(svd1$d, ylab = "Singular value", col = 2, main = "Singular Values")
dev.off()

# Plot percent variance explained by each component
png("percent_variance_explained.png", width = 1600, height = 1600, res = 250)
plot(svd1$d^2 / sum(svd1$d^2), ylab = "Percent Variance Explained", col = 2, main = "Percent Variance Explained")
dev.off()

# Plot first two principal components
png("pc1_pc2.png", width = 1600, height = 1600, res = 250)
par(mfrow = c(1, 2))
plot(svd1$v[, 1], col = 2, ylab = "1st PC", main = "First Principal Component")
plot(svd1$v[, 2], col = 2, ylab = "2nd PC", main = "Second Principal Component")
dev.off()

# PC1 vs PC2 scatter plot
png("pc1_vs_pc2.png", width = 1600, height = 1600, res = 250)
plot(svd1$v[, 1], svd1$v[, 2], col = 2, ylab = "2nd PC", xlab = "1st PC", main = "PC1 vs PC2")
dev.off()

# Color points by study
png("pc1_vs_pc2_by_study.png", width = 1600, height = 1600, res = 250)
plot(svd1$v[, 1], svd1$v[, 2], ylab = "2nd PC", xlab = "1st PC", col = as.numeric(pdata$study), main = "PC1 vs PC2 (Colored by Study)")
dev.off()

# Compare PC distributions between studies
png("pc1_distribution_by_study.png", width = 1600, height = 1600, res = 250)
boxplot(svd1$v[, 1] ~ pdata$study, border = c(1, 2), main = "PC1 Distribution by Study")
points(svd1$v[, 1] ~ jitter(as.numeric(pdata$study)), col = as.numeric(pdata$study))
dev.off()

# Compare with standard PCA
pc1 = prcomp(edata)
png("svd_vs_pca.png", width = 1600, height = 1600, res = 250)
plot(pc1$rotation[, 1], svd1$v[, 1], col = 2, xlab = "PCA Rotation", ylab = "SVD V", main = "SVD vs PCA")
dev.off()

# Column-centered SVD
edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)
png("svd_column_centered.png", width = 1600, height = 1600, res = 250)
plot(pc1$rotation[, 1], svd2$v[, 1], col = 2, xlab = "PCA Rotation", ylab = "SVD V (Column Centered)", main = "SVD Column Centered")
dev.off()

# Create and analyze outlier
edata_outlier = edata_centered
edata_outlier[1, ] = edata_centered[1, ] * 10000
svd3 = svd(edata_outlier)

# Compare SVD with/without outlier
png("svd_with_without_outlier.png", width = 1600, height = 1600, res = 250)
par(mfrow = c(1, 2))
plot(svd1$v[, 1], col = 1, main = "Without Outlier")
plot(svd3$v[, 1], col = 2, main = "With Outlier")
dev.off()

# Plot correlation between SVD and outlier
png("svd_outlier_correlation.png", width = 1600, height = 1600, res = 250)
plot(svd3$v[, 1], edata_outlier[1, ], col = 4, xlab = "SVD V", ylab = "Outlier Expression", main = "SVD vs Outlier")
dev.off()

