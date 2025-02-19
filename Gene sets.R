# Load required libraries
library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "goseq", "DESeq2"))

# Download the data
# Check supported genomes and gene IDs by goseq
head(supportedGenomes())
head(supportedGeneIDs())

# An example of a goseq analysis
# Load the data
temp_data = read.table(system.file("extdata", "Li_sum.txt", package="goseq"), 
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
expr = temp_data[,-1]
rownames(expr) = temp_data[,1]
expr = expr[rowMeans(expr) > 5,]  # Filter out low expression genes
grp = factor(rep(c("Control", "Treated"), times=c(4,3)))
pdata = data.frame(grp)

# Perform a differential expression analysis with DESeq2
de = DESeqDataSetFromMatrix(expr, pdata, ~grp)
de_fit = DESeq(de)
de_results = results(de_fit)

# Get the differentially expressed genes after FDR correction
genes = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes)
names(genes) = rownames(expr)
genes = genes[not_na]

# Pick the right genome
head(supportedGenomes(), n=12)[,1:4]

# Set up a weighting function for all genes in the genome
pwf = nullp(genes, "hg19", "ensGene")
head(pwf)

# Perform the enrichment analysis parametrically
GO.wall = goseq(pwf, "hg19", "ensGene")
head(GO.wall)

# Limiting to a single category of interest
GO.MF = goseq(pwf, "hg19", "ensGene", test.cats=c("GO:MF"))
head(GO.MF)
