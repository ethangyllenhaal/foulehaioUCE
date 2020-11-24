##########
# By: Ethan Gyllenhaal
# Updated 16Nov2020
#
# R script used for generating PCAs for Foulehaio carunculatus.
# First initializes PCA colors, then makes a PCA for each one.
# Then has one section for each relevant PCA.
# Per-section process described in the first section.
# Note that PDFs were made manually, not output to a file.
########

# Packages! adegenet and ade4 are for the PCAs, vcfR is file conversion, scales is colors
library("adegenet")
library("ade4")
library("vcfR")
library("scales")
setwd('/path/to/directory')

# sets colors
upolu = rgb(0, 220, 180, 255, maxColorValue = 255)
sav = rgb(0, 200, 200, 255, maxColorValue = 255)
tut = rgb(0, 150, 200, 255, maxColorValue = 255)
tau = rgb(0, 40, 200, 255, maxColorValue = 255)
olo = rgb(0, 150, 255, 255, maxColorValue = 255)
ofu = rgb(0, 100, 255, 255, maxColorValue = 255)
eueki = rgb(20, 255, 100, 255, maxColorValue = 255)
tofua = rgb(40, 180, 90, 255, maxColorValue = 255)
eua = rgb(60, 150, 60, 255, maxColorValue = 255)
alofi = rgb(190, 210, 170, 255, maxColorValue = 255)
futuna = rgb(120, 170, 120, 255, maxColorValue = 255)
ono = rgb(120, 30, 120, 255, maxColorValue = 255)
matuku = rgb(80, 10, 80, 255, maxColorValue = 255)
fulaga = rgb(170, 120, 180, 255, maxColorValue = 255)
ogea = rgb(210, 140, 220, 255, maxColorValue = 255)

### Samoa and Tonga (explanation of each section here)
# read in VCF
ST_vcf <- read.vcfR("carunculatus_samoa-tonga_PCA-input.vcf")
# converts VCF to genlight format
ST_gl <- vcfR2genlight(ST_vcf, n.cores = 1)
# assigns populations to samples, in VCF order
pop(ST_gl) <- c("'Eua", "Savai'i", "Savai'i", "Upolu", "Upolu", "Eueki", "Upolu", "Tau", "Tutuila", "Tutuila", "Tutuila", "Tutuila", "Olosega", "Eueki", "Eueki", "Eueki", "Eueki", "Tofua", "Tofua", "Olosega", "Savai'i", "Savai'i", "Savai'i", "Savai'i", "Upolu", "Upolu", "Olosega", "Ofu", "Ofu", "Tutuila", "Tau", "Tau", "Ofu", "Ofu", "Upolu")
# runs PCA, retaining "nf" components
pca_ST <- glPca(ST_gl, nf=20)
# plots PCA with colors
s.class(pca_ST$scores, pop(ST_gl), col=c(eua, eueki, ofu, olo, sav, tau, tofua, tut, upolu), clab=1, cstar=0, cell=2.5)
# plots with PC3 and 4
s.class(pca_ST$scores[,3:4], pop(ST_gl), col=c(eua, eueki, ofu, olo, sav, tau, tofua, tut, upolu), clab=1, cstar=0, cell=2.5)
# makes an optional barplot of eigenvalues to evaluate effect of each PC
barplot(pca_ST$eig/sum(pca_ST$eig), main="eigenvalues", col=heat.colors(length(pca_ST$eig)))


# Tonga
tonga_vcf <- read.vcfR("carunculatus_tonga_PCA-input.vcf")
tonga_gl <- vcfR2genlight(tonga_vcf, n.cores = 1)
pop(tonga_gl) <- c("'Eua", "Eueki", "Eueki", "Eueki", "Eueki", "Eueki", "Tofua", "Tofua")
pca_tonga <- glPca(tonga_gl, nf=20)
s.class(pca_tonga$scores, pop(tonga_gl), col=c(eua, eueki, tofua), clab=1, cstar=0, cell=2.5)
s.class(pca_tonga$scores[,2:3], pop(tonga_gl), col=c("blue","green","black","orange","red","purple", "black", "brown", "gray"), clab=1, cstar=0, cell=2.5)
barplot(pca_tonga$eig/sum(pca_tonga$eig), main="eigenvalues", col=heat.colors(length(pca_tonga$eig)))

# Samoa
samoa_vcf <- read.vcfR("carunculatus_samoa_PCA-input.vcf")
samoa_gl <- vcfR2genlight(samoa_vcf, n.cores = 1)
pop(samoa_gl) <- c("Savai'i", "Savai'i", "Upolu", "Upolu", "Upolu", "Tau", "Tutuila", "Tutuila", "Tutuila", "Tutuila", "Olosega", "Olosega", "Savai'i", "Savai'i", "Savai'i", "Savai'i", "Upolu", "Upolu", "Olosega", "Ofu", "Ofu", "Tutuila", "Tau", "Tau", "Ofu", "Ofu", "Upolu")
pca_samoa <- glPca(samoa_gl, nf=20)
s.class(pca_samoa$scores, pop(samoa_gl), col=c(ofu, olo, sav, tau, tut, upolu), clab=1, cstar=0, cell=2.5)
s.class(pca_samoa$scores[,3:4], pop(samoa_gl), col=c(ofu, olo, sav, tau, tut, upolu), clab=1, cstar=0, cell=2.5)
barplot(pca_samoa$eig/sum(pca_samoa$eig), main="eigenvalues", col=heat.colors(length(pca_samoa$eig)))

# Samoa WC
samoaWC_vcf <- read.vcfR("carunculatus_samoa-WandC_PCA-input.vcf")
samoaWC_gl <- vcfR2genlight(samoaWC_vcf, n.cores = 1)
pop(samoaWC_gl) <- c("Savai'i", "Savai'i", "Upolu", "Upolu", "Upolu", "Tutuila", "Tutuila", "Tutuila", "Tutuila", "Savai'i", "Savai'i", "Savai'i", "Savai'i", "Upolu", "Upolu", "Tutuila", "Upolu")
pca_samoaWC <- glPca(samoaWC_gl, nf=20)
s.class(pca_samoaWC$scores, pop(samoaWC_gl), col=c(sav, tut, upolu), clab=1, cstar=0, cell=2.5)
barplot(pca_samoaWC$eig/sum(pca_samoaWC$eig), main="eigenvalues", col=heat.colors(length(pca_samoaWC$eig)))

# Samoa E
samoaE_vcf <- read.vcfR("carunculatus_samoa-E_PCA-input.vcf")
samoaE_gl <- vcfR2genlight(samoaE_vcf, n.cores = 1)
pop(samoaE_gl) <- c("Tau", "Olosega", "Olosega", "Olosega", "Ofu", "Ofu", "Tau", "Tau", "Ofu", "Ofu")
pca_samoaE <- glPca(samoaE_gl, nf=20)
s.class(pca_samoaE$scores[,1:2], pop(samoaE_gl), col=c(ofu, olo, tau), clab=1, cstar=0, cell=2.5)
barplot(pca_samoaE$eig/sum(pca_samoaE$eig), main="eigenvalues", col=heat.colors(length(pca_samoaE$eig)))

# Alofi/Futuna
AF_vcf <- read.vcfR("carunculatus_futuna_PCA-input.vcf")
AF_gl <- vcfR2genlight(AF_vcf, n.cores = 1)
pop(AF_gl) <- c("Alofi", "Alofi", "Futuna", "Futuna", "Futuna", "Futuna", "Futuna", "Alofi", "Alofi", "Alofi")
pca_AF <- glPca(AF_gl, nf=20)
s.class(pca_AF$scores, pop(AF_gl), col=c(alofi, futuna), clab=1, cstar=0, cell=2.5)
s.class(pca_AF$scores[,2:3], pop(AF_gl), col=c(alofi, futuna), clab=1, cstar=0, cell=2.5)
barplot(pca_AF$eig/sum(pca_AF$eig), main="eigenvalues", col=heat.colors(length(pca_AF$eig)))

# Lau
lau_vcf <- read.vcfR("carunculatus_lau_PCA-input.vcf")
lau_gl <- vcfR2genlight(lau_vcf, n.cores = 1)
pop(lau_gl) <- c("Ono-i-Lau", "Ono-i-Lau", "Ono-i-Lau", "Ono-i-Lau", "Ogea Levu", "Ogea Levu", "Ogea Levu", "Fulaga", "Fulaga", "Fulaga", "Fulaga", "Matuku")
pca_lau <- glPca(lau_gl, nf=20)
s.class(pca_lau$scores, pop(lau_gl), col=c(fulaga, matuku, ogea, ono), clab=1, cstar=0, cell=2.5)
s.class(pca_lau$scores[,2:3], pop(lau_gl), col=c(fulaga, matuku, ogea, ono), clab=1, cstar=0, cell=2.5)
barplot(pca_lau$eig/sum(pca_lau$eig), main="eigenvalues", col=heat.colors(length(pca_lau$eig)))
