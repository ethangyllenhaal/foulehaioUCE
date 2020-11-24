###################
#
# By: Ethan Gyllenhaal
# Last updated: 18Nov2020
#
# Script for running ABBA/BABA with bootstrapping for both D and f statistics.
# General framework inspired by and using scripts by Simon Martin: https://github.com/simonhmartin/genomics_general
# First batch of functions are simply calculating and printing the test statistics.
# Second batch (calc_*) of functions are drivers for bootstrapping and making the density curves.
# Third batch has "bootstrap_*" and "bootrep_*" functions, for generating bootstrapped replicates and summarizing them.
# Final output is the test stat, a Z value, a bootstrapped p value (# of reps greater than 0), and a density curve of the statistic.
#
##################

# Function to calculate D statistic given three columns of a derived frequency table.
D.stat <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
}

# Prints out the values of ABBA and BABA used for calculation, then D statistic.
D.stat_print <- function(p1, p2, p3) {
  ABBA <- (1 - p1) * p2 * p3
  BABA <- p1 * (1 - p2) * p3
  print(paste("ABBA:",sum(ABBA),"BABA:", sum(BABA)))
  print(paste("D statistic:", (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))))
}

# Function to calculate f given two unique columns of a derived frequency table and two "P3" pops.
# Populations p3a and p3b can be subpopulations, arbitrarily defined pops, or the same pop (as used below).
f.stat <- function(p1, p2, p3a, p3b) {
  ABBA_numerator <- (1 - p1) * p2 * p3a
  BABA_numerator <- p1 * (1 - p2) * p3a
  
  ABBA_denominator <- (1 - p1) * p3b * p3a
  BABA_denominator <- p1 * (1 - p3b) * p3a
  
  (sum(ABBA_numerator) - sum(BABA_numerator)) /
    (sum(ABBA_denominator) - sum(BABA_denominator))
}

# Driver function for D stat here, calls other functions and gives them proper input.
# "in_name" is frequency table, P1-3 are strings of name in frequency table, and bootnum is # of bootstraps to perform.
# Outputs ABBA/BABA counts, the D statistic, then bootstrapped p value (i.e. count(D<=0)/bootreps), and finally Z score.
calc_D_stat <- function(in_name, P1, P2, P3, bootnum){
  D <- D.stat(in_name[,P1], in_name[,P2], in_name[,P3])
  D.stat_print(in_name[,P1], in_name[,P2], in_name[,P3])
  block_indices <- get_block_indices(block_size=1e6, positions=in_name$position, chromosomes=in_name$scaffold)
  n_blocks <- length(block_indices)
  D_boot <- bootstrap_D(block_indices, in_name[,P1], in_name[,P2], in_name[,P3], bootnum)
  D_Z <- D / sd(D_boot)
  print(paste("The Z score is", D_Z))
  D_boot
}

# Driver function for f-hom here, calls other functions and gives them proper input.
# "in_name" is frequency table, P1-3 are strings of name in frequency table, and bootnum is # of bootstraps to perform.
# Outputs fhom, then bootstrapped p value (i.e. count(D<=0)/bootreps), and 95% Z confidence interval of fhom.
calc_f_hom <- function(in_name, P1, P2, P3, bootnum){
  f <- f.stat(in_name[,P1], in_name[,P2], in_name[,P3], in_name[,P3])
  print(paste("F hom:",f))
  # get block indices from jacknife.R in https://github.com/simonhmartin/genomics_general
  block_indices <- get_block_indices(block_size=1e6, positions=in_name$position, chromosomes=in_name$scaffold)
  n_blocks <- length(block_indices)
  f_boot <- bootstrap_f(block_indices, in_name[,P1], in_name[,P2], in_name[,P3], in_name[,P3], bootnum)
  f_sd <- sd(f_boot)
  f_CI_lower <- f - 1.96*f_sd
  f_CI_upper <- f + 1.96*f_sd
  print(paste("95% confidence interval of f =", round(f_CI_lower,4), round(f_CI_upper,4)))
  f_boot
}

###################
## Bootstrapping ##
###################

# Main function for doing all bootstaps.
# Input is block indices, P1-3, and number of reps.
bootstrap_D <- function(index, P1, P2, P3, bootnum){
  arr = c()
  # runs bootstraps
  for (i in 1:bootnum){
    arr <- c(arr, bootrep_D(index, P1, P2, P3))
  }
  sum_neg = 0
  # counts number of reps with D <=0
  for (j in arr){
    if (j <= 0){
      sum_neg <- sum_neg + 1
    }
  }
  print(paste(c("Proportion of values below zero:", sum_neg/bootnum), collapse = " "))
  d <- density(arr)
  # makes density plot, you likely want to change x and y limits as appropriate
  x_label <- paste("D = ", round(D.stat(P1, P2, P3), digits=3), " p = ", round(sum_neg/bootnum, digits=3))
  # Commented out plot function is limits used in paper, can change as needed for clarity and consistency
  #plot(d, main = "", ylab = "", xlab = x_label, xlim = c(-0.7, 0.7), ylim=c(0,3))
  plot(d, main = "", ylab = "", xlab = x_label)
  polygon(d, col="gray")
  abline(v=0)
  arr
}

# Function for a single bootrep of D.
# Input is block indices and P1-3.
bootrep_D <- function(index, P1, P2, P3){
  p1_boot = c()
  p2_boot = c()
  p3_boot = c()
  len = length(index)
  counter = 0
  list = list()
  # goes through and adds relevant sites until array full (to account for freq table w/ more than 3 taxa)
  while(counter != len){
    temp <- sample(1:length(P1), 1, replace=TRUE)
    # only adds sites with frequency of P1, P2, and/or P3 present
    if(P1[temp]+P2[temp]+P3[temp]>0){
      list <- c(list, temp)
      counter = counter + 1     
    }
  }
  count = 0
  for (i in list){
    p1_boot <- c(p1_boot, P1[i])
    p2_boot <- c(p2_boot, P2[i])
    p3_boot <- c(p3_boot, P3[i])
    
  }
  D.stat(p1_boot, p2_boot, p3_boot)
}

# Main function for doing all bootstaps of fhom.
# Input is block indices, P1-3, and number of reps.
bootstrap_f <- function(index, P1, P2, P3a, P3b, bootnum){
  arr = c()
  # runs bootstraps
  for (i in 1:bootnum){
    arr <- c(arr, bootrep_f(index, P1, P2, P3a, P3b))
  }
  sum_neg = 0
  # counts number of reps with D <=0
  for (j in arr){
    if (j <= 0){
      sum_neg <- sum_neg + 1
    }
  }
  print(paste(c("Proportion of values below zero:", sum_neg/bootnum), collapse = " "))
  # makes density plot, you likely want to change x and y limits as appropriate
  f <- density(arr)
  x_label <- paste("f = ", round(f.stat(P1, P2, P3a, P3b), digits=3), " p = ", round(sum_neg/bootnum, digits=3))
  # Commented out plot function is limits used in paper, can change as needed for clarity and consistency
  #plot(f, main = "", ylab = "", xlab = x_label, xlim = c(-0.1, 0.1), ylim=c(0,25))
  plot(f, main = "", ylab = "", xlab = x_label)
  polygon(f, col="gray")
  abline(v=0)
  arr
}

# Function for a single bootrep of fhom.
# Input is block indices and P1-3b.
bootrep_f <- function(index, P1, P2, P3a, P3b){
  p1_boot = c()
  p2_boot = c()
  p3a_boot = c()
  p3b_boot = c()
  len = length(index)
  list = list()
  counter = 0
  # only adds sites with frequency of P1, P2, and/or P3 present
  while(counter != len){
    temp <- sample(1:length(P1), 1, replace=TRUE)
    if(P1[temp]+P2[temp]+P3a[temp]>0){
      list <- c(list, temp)
      counter = counter + 1     
    }
  }
  for (i in list){
    p1_boot <- c(p1_boot, P1[i])
    p2_boot <- c(p2_boot, P2[i])
    p3a_boot <- c(p3a_boot, P3a[i])
    p3b_boot <- c(p3b_boot, P3b[i])
  }
  f.stat(p1_boot, p2_boot, p3a_boot, p3b_boot)
}
