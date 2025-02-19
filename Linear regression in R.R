# Load required libraries
library(devtools)
library(Biobase)
library(broom)

# Install necessary packages if not already installed
# Uncomment these lines to install packages
# install.packages(c("devtools","broom"))
# source("http://www.bioconductor.org/biocLite.R")
# biocLite(c("Biobase"))

# Download and load the data
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)  # Phenotypic data
edata = as.data.frame(exprs(bm))  # Expression data
fdata = fData(bm)  # Feature data

# Fit a simple linear regression for the first gene against age
edata = as.matrix(edata)
lm1 = lm(edata[1,] ~ pdata$age)
tidy(lm1)  # Display summary of linear model

# Visual diagnostics for the linear regression
plot(pdata$age, edata[1,], col=1)
abline(lm1$coeff[1], lm1$coeff[2], col=2, lwd=3)

# Visualize the difference between genders for the first gene
boxplot(edata[1,] ~ pdata$gender)
points(edata[1,] ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender))

# Create dummy variables for gender
dummy_m = pdata$gender == "M"
dummy_f = pdata$gender == "F"

# Fit linear model with gender
lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)  # Display summary of linear model with gender

# Check model matrix for gender
mod2 = model.matrix(~pdata$gender)
mod2

# Explore categorical variables with multiple levels
table(pdata$tissue.type)
# Check for specific tissue types
pdata$tissue.type == "adipose"
pdata$tissue.type == "adrenal"

# Linear model with tissue type
tidy(lm(edata[1,] ~ pdata$tissue.type))

# Adjusting for both age and gender
lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)

# Interaction effects between age and gender
lm4 = lm(edata[1,] ~ pdata$age * pdata$gender)
tidy(lm4)

# Example of outlier impact on regression
lm4 = lm(edata[6,] ~ pdata$age)
plot(pdata$age, edata[6,], col=2)
abline(lm4, col=1, lwd=3)

# Another example with an obvious outlier
index = 1:19
lm5 = lm(edata[6,] ~ index)
plot(index, edata[6,], col=2)
abline(lm5, col=1, lwd=3)

# Remove outlier and compare
lm6 = lm(edata[6,-19] ~ index[-19])
abline(lm6, col=3, lwd=3)
legend(5, 1000, c("With outlier", "Without outlier"), col=c(1,3), lwd=3)

# Check residuals distribution
par(mfrow=c(1,2))
hist(lm6$residuals, col=2)
hist(lm5$residuals, col=3)

# Data transformation and residuals
gene1 = log2(edata[1,] + 1)  # Log transform data
lm7 = lm(gene1 ~ index)
hist(lm7$residuals, col=4)

# Example of co-linearity problem
lm8 = lm(gene1 ~ pdata$tissue.type + pdata$age)
tidy(lm8)  # Note: This might show NaN due to co-linearity

# Check for patterns in residuals
colramp = colorRampPalette(1:4)(17)
lm9 = lm(edata[2,] ~ pdata$age)
plot(lm9$residuals, col=colramp[as.numeric(pdata$tissue.type)])

