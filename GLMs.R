# Load required libraries
library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)

# Install necessary packages if not already installed
# Uncomment these lines to install packages
# install.packages(c("devtools", "broom", "MASS"))
# source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("Biobase", "snpStats", "DESeq2"))

# Logistic regression

# Load the data
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]

# Calculate the PCs
xxmat <- xxt(sub.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]

# A single logistic regression
snpdata = sub.10@.Data
status = subject.support$cc
snp1 = as.numeric(snpdata[,1])
snp1[snp1==0] = NA
glm1 = glm(status ~ snp1, family="binomial")
tidy(glm1)

# Dominant model
snp1_dom = (snp1 == 1)
glm1_dom = glm(status ~ snp1_dom, family="binomial")
tidy(glm1_dom)

# Adjust for principal components
glm2 = glm(status ~ snp1 + pcs[,1:5], family="binomial")
tidy(glm2)

# Fit many GLMs at once
glm_all = snp.rhs.tests(status ~ 1, snp.data=sub.10)
slotNames(glm_all)
qq.chisq(chi.squared(glm_all), df=1)

# Adjust for PCs in many GLMs
glm_all_adj = snp.rhs.tests(status ~ pcs, snp.data=sub.10)
qq.chisq(chi.squared(glm_all_adj), df=1)

# Poisson/negative binomial regression

# Download the data
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata = pData(bot)
edata = as.matrix(exprs(bot))
fdata = fData(bot)

# Transform the data
edata = edata[rowMeans(edata) > 10, ]

# A single Poisson regression
glm3 = glm(edata[1, ] ~ pdata$strain, family="poisson")
tidy(glm3)

# A single negative binomial regression
glm.nb1 = glm.nb(edata[1, ] ~ pdata$strain)
tidy(glm.nb1)

# Multiple negative binomial regressions with DESeq2
de = DESeqDataSetFromMatrix(edata, pdata, ~strain)
glm_all_nb = DESeq(de)
result_nb = results(glm_all_nb)
hist(result_nb$stat)
