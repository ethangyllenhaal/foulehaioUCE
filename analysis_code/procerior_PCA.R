##########
# By: Ethan Gyllenhaal
# Updated 16Nov2020
#
# R script used for generating PCAs for Foulehaio procerior.
# First initializes PCA colors, then makes a PCA for the species.
# Note that PDF was made manually, not output to a file.
########

# Packages! adegenet and ade4 are for the PCAs, vcfR is file conversion.
library("adegenet")
library("ade4")
library("vcfR")
setwd('/path/to/directory')

# sets colors
vit1 = rgb(200, 0, 0, 255, maxColorValue = 255)
vit2 = rgb(255, 50, 0, 255, maxColorValue = 255)

# read in VCF
car_vcf <- read.vcfR("procerior_PCA-input.vcf")
# converts VCF to genlight format
car.genlight <- vcfR2genlight(car_vcf, n.cores=1)
# assigns populations to samples, in VCF order
pop(car.genlight) <- c("viti", "viti", "viti", "viti", "ovalau", "ovalau", "ovalau", "viti", "viti")
# runs PCA, retaining "nf" components
pca1 <- glPca(car.genlight, nf=20)
# plots PCA with colors
s.class(pca1$scores, pop(car.genlight), col=c(vit1, vit2, "blue", "black"), clab=1, cell=2.5)
# makes an optional barplot of eigenvalues to evaluate effect of each PC
barplot(pca1$eig/sum(pca1$eig), main="eigenvalues", col=heat.colors(length(pca1$eig)))

