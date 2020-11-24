###################
#
# By: Ethan Gyllenhaal
# Last updated: 21Nov2020
#
# Script for running sNMF for F. carunculatus, all populationjs together.
# First sets up sNMF driver function, then runs sNMF, and finally plots for different K values.
#
##################

# Load in additional tools for sNMF
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
# Load main library, LEA
library(LEA)
setwd('path/to/directory/')

# Plotting function used to plot sNMF output.
# Input is sNMF object, the k value, and optionally an array of colors (has default of 10).
plot_sNMF <- function(input, k_val, colors = c("green", "red", "yellow", "blue", "orange", "purple", "tan", "gray", "dark green", "pink")){
  # picks best run best on cross entropy
  best_run <- which.min(cross.entropy(input, K = k_val))
  # makes q matrix of ancestry coeffs
  q_matrix <- Q(input, K = k_val, run = best_run)
  # plots the output, space makes blank between indivs
  barplot(t(q_matrix), col = colors, border = NA, space = 0.25, xlab = "Individuals", ylab = "Admixture coefficients", horiz=FALSE)
}


# convert .str file to a .geno file, with first column being sample name and second being population #
# If no population column (they aren't really used anyways), set extra.col to 0
struct2geno(file = "carunc_100_rand_UG_noSingle.str", 
            TESS = FALSE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 1, output = "D:/Documents/Projects/Foulehaio/carunc_redo/sNMF_carunc_redo.geno")

# Runs sNMF for all values of K (K=1:12 in my case) and 10 repetitions. The entropy = T is how you estimate best K.
# Note that alpha is the "cost" of introducing admixture, and higher values = less admixture
carunc_redo = snmf("D:/Documents/Projects/Foulehaio/carunc_redo/sNMF_carunc_redo.geno", ploidy=2, 
                K = 1:12, alpha = 10, project = "new", entropy = T, repetitions = 10)

# plots values of K, choose K as per structure 
plot(carunc_redo, cex = 1.2, col = "lightblue", pch = 19)

# plotting a few values of k
# 1-4 Onoilau, 5-11 Lau, 12 Matuku, 13-20 Alofi/Futuna, 21-28 Tonga, 29-36 Samoa
# 5-7 and 8-11 diff islands in Lau, 21-25 and 26-28 different parts of tonga, 29-32 and 33-36 different parts of samoa
plot_sNMF(carunc_redo, 3)
plot_sNMF(carunc_redo, 4)
plot_sNMF(carunc_redo, 5)
plot_sNMF(carunc_redo, 6)
plot_sNMF(carunc_redo, 7)
plot_sNMF(carunc_redo, 8)
plot_sNMF(carunc_redo, 9)
plot_sNMF(carunc_redo, 10)
plot_sNMF(carunc_redo, 11)
plot_sNMF(carunc_redo, 12)


