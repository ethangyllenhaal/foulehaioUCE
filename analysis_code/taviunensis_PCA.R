##########
# By: Ethan Gyllenhaal
# Updated 16Nov2020
#
# R script used for generating PCAs for Foulehaio taviunensis.
# First initializes PCA colors, then makes a PCA for the species.
# Note that PDF was made manually, not output to a file.
########

# Packages! adegenet and ade4 are for the PCAs, vcfR is file conversion.
library("adegenet")
library("ade4")
library("vcfR")
setwd('/path/to/directory')

# sets colors
tav = rgb(200, 0, 0, 255, maxColorValue = 255)
van = rgb(255, 50, 0, 255, maxColorValue = 255)
kio = rgb(255, 150, 0, 255, maxColorValue = 255)

# read in VCF
car_vcf <- read.vcfR("taviunensis_PCA-input.vcf")
# converts VCF to genlight format
car.genlight <- vcfR2genlight(car_vcf, n.cores=1)
# assigns populations to samples, in VCF order
pop(car.genlight) <- c("Vanua Levu", "Vanua Levu", "Vanua Levu", "Kioa", "Vanua Levu", "Vanua Levu", "Vanua Levu", "Vanua Levu", "Taveuni", "Taveuni", "Taveuni")
# runs PCA, retaining "nf" components
pca1 <- glPca(car.genlight, nf=20)
# plots PCA with colors
s.class(pca1$scores, pop(car.genlight), col=c(kio,tav,van), clab=1, cell=2.5)
# plots with PC2 and 3
s.class(pca1$scores[,2:3], pop(car.genlight), col=c(kio,tav,van), clab=1, cell=2.5)
# makes an optional barplot of eigenvalues to evaluate effect of each PC
barplot(pca1$eig/sum(pca1$eig), main="eigenvalues", col=heat.colors(length(pca1$eig)))

