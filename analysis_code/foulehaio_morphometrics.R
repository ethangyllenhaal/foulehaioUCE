###################
#
# By: Ethan Gyllenhaal
# Last updated: 23Nov2020
#
# Script where morphometric analyses were run from, first RAXML then SNAPP topologies.
# Each section first loads and prunes tree, prepares the dataset, tests for phylogenetic signal, and runs PGLS.
# Finishes making ContMaps (trait mapped to tree) for both topologies.
#
##################

# loads in a lot of packages, pretty sure all of these are used in final version of script?
library(ape)
library(data.table)
library(phytools)
library(phylotools)
library(geiger)
library(picante)
library(nlme)
library(ggplot2)
library(caper)
library(wiqid)
# set path to working directory
setwd('path/to/directory')



###########
## RAxML ##
###########

# loads in full RAxML tree
full_tree <- read.tree("foulehaio_raxml_tree.tre")

# makes list of tips to keep
keepers <- c("foulehaio_procerior_22540", "foulehaio_taviunensis_24322", "foulehaio_taviunensis_waho_fjos_85_val", "foulehaio_carunculatus_26425", "foulehaio_carunculatus_26309", "foulehaio_carunculatus_26352", "foulehaio_carunculatus_26385", "foulehaio_carunculatus_waho_2012_162_al", "foulehaio_carunculatus_rbu_2012_192", "foulehaio_carunculatus_rbu_2012_356", "foulehaio_carunculatus_jos_sa_142_waho", "foulehaio_carunculatus_rbu_2012_053", "foulehaio_carunculatus_waho_311_tau")
# reads in data .csv for PGLS
dataset <- read.csv("combined_measures_adjusted_full.csv", header=TRUE, sep=",", row.names = 1)
# trims tree
trimmed_tree <- keep.tip(full_tree, keepers)
# plots tree to check
plot(root(trimmed_tree, "foulehaio_procerior_22540"))

# adjusts labels, makes tree ultrametric
final_tree <- sub.taxa.label(trimmed_tree, as.data.frame(cbind(dataset["Sample"][,1], row.names(dataset))))
final_tree <- root(final_tree, "Procerior", resolve.root = TRUE)
final_tree <- reroot(final_tree, 1, position=0.1 * final_tree$edge.length[which(final_tree$edge[,2] == 1)])
final_tree <- force.ultrametric(final_tree, method=("extend"))

# takes out male and female data
males <- dataset["Male"]
males <- as.matrix(males)[,1]
females <- dataset["Female"]
females <- as.matrix(females)[,1]

# runs both tests of plohygenetic signal
phylosig(final_tree, males, method="K", test=TRUE, nsim=1000)
phylosig(final_tree, males, method="lambda", test=TRUE)


# sets PGLS tree to island-specific carunculatus data only
pgls_tree <- drop.tip(final_tree, c("Procerior","Taveuni", "Vanua", "Tonga"))
# generates correlation matrix
brownian <- corBrownian(1, pgls_tree)

# initializes all the different models
model_int <- update(gls(Male ~ Lat * log(Size), data=dataset[4:12,], correlation = brownian), . ~ ., method = "ML")
model_add <- update(gls(Male ~ Lat + log(Size), data=dataset[4:12,], correlation = brownian), . ~ ., method = "ML")
model_size <- update(gls(Male ~ log(Size), data=dataset[4:12,], correlation = brownian), . ~ ., method = "ML")
model_lat <- update(gls(Male ~ Lat, data=dataset[4:12,], correlation = brownian), . ~ ., method = "ML")
model_null <- update(gls(Male ~ 1, data=dataset[4:12,], correlation = brownian), . ~ ., method = "ML")

# gets AIC and p-value for each model
# (I didn't set this up to populate a table automatically, kind of done the "dumb" way)
summary(model_int)$AIC
summary(model_int)$t
summary(model_add)$AIC
summary(model_add)$t
summary(model_size)$AIC
summary(model_size)$t
summary(model_lat)$AIC
summary(model_lat)$t
summary(model_null)$AIC
summary(model_null)$t

###########
## SNAPP ##
###########

# load in SNAPP tree and dataset without Matuku
dataset_snapp <- read.csv("combined_measures_adjusted_noMat.csv", header=TRUE, sep=",", row.names = 1)
snapp_tree <- read.nexus("foulehaio_island_snapp_tree.tre")

# males male and female datasets
snapp_m <- dataset_snapp["Male"]
snapp_m <- as.matrix(snapp_m)[,1]
snapp_f <- dataset_snapp["Female"]
snapp_f <- as.matrix(snapp_f)[,1]

# tests for phylogenetic signal
phylosig(snapp_tree, snapp_m, method="K", test=TRUE, nsim=1000)
phylosig(snapp_tree, snapp_m, method="lambda", test=TRUE)

# sets PGLS tree to island-specific carunculatus data only
pgls_snapp <- drop.tip(snapp_tree, c("Vanua", "Tonga"))
# makes correlation matrix
brownian_snapp <- corBrownian(1, pgls_snapp)


# initializes all the different models
snapp_int <- update(gls(Male ~ Lat * log(Size), data=dataset_snapp[2:9,], correlation = brownian_snapp), . ~ ., method = "ML")
snapp_add <- update(gls(Male ~ Lat + log(Size), data=dataset_snapp[2:9,], correlation = brownian_snapp), . ~ ., method = "ML")
snapp_size <- update(gls(Male ~ log(Size), data=dataset_snapp[2:9,], correlation = brownian_snapp), . ~ ., method = "ML")
snapp_lat <- update(gls(Male ~ Lat, data=dataset_snapp[2:9,], correlation = brownian_snapp), . ~ ., method = "ML")
snapp_null <- update(gls(Male ~ 1, data=dataset_snapp[2:9,], correlation = brownian_snapp), . ~ ., method = "ML")

# gets AIC and p-values for each model
# (I didn't set this up to populate a table automatically, kind of done the "dumb" way)
summary(snapp_int)$AIC
summary(snapp_int)$t
summary(snapp_add)$AIC
summary(snapp_add)$t
summary(snapp_size)$AIC
summary(snapp_size)$t
summary(snapp_lat)$AIC
summary(snapp_lat)$t
summary(snapp_null)$AIC
summary(snapp_null)$t


#############
## ContMap ##
#############

# rotates trees to make them match
holder_tree <- final_tree
final_tree <- rotateNodes(holder_tree, c(16, 15, 18, 21))
snapp_tree <- rotateNodes(snapp_tree, c(13, 14, 16, 17))

# plots to check concordance
par(mfrow=c(1,2))
plotTree(final_tree, node.numbers=T)
plotTree(snapp_tree, node.numbers=T)

# plots each ContMap by female and tree type
par(mfrow=c(2,2))
cont_rax_male <- contMap(final_tree, males, fsize=c(1,0.6), plot=FALSE)
plot(setMap(cont_rax_male, colors=c("blue", "cyan", "lawngreen", "yellow", "red")), outline=FALSE)
title(main = "Male Wing Length - RAxML", cex.main=0.7)
cont_rax_female <- contMap(final_tree, females, fsize=c(1,0.6), plot=FALSE)
plot(setMap(cont_rax_female, colors=c("blue", "cyan", "lawngreen", "yellow", "red")), outline=FALSE)
title(main = "Female Wing Length - RAxML", cex.main=0.7)
cont_snapp_male <- contMap(snapp_tree, snapp_m, fsize=c(1,0.6), plot=FALSE)
plot(setMap(cont_snapp_male, colors=c("blue", "cyan", "lawngreen", "yellow", "red")), outline=FALSE)
title(main = "Male Wing Length - SNAPP", cex.main=0.7)
cont_snapp_female <- contMap(snapp_tree, snapp_f, fsize=c(1,0.6), plot=FALSE)
plot(setMap(cont_snapp_female, colors=c("blue", "cyan", "lawngreen", "yellow", "red")), outline=FALSE)
title(main = "Female Wing Length - SNAPP", cex.main=0.7)