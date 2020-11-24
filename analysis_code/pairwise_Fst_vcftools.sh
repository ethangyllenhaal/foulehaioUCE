#!/bin/bash

# Made by Ethan F. Gyllenhaal (egyllenhaal@unm.edu)
# Last updated: 13Nov2020
#
# This is a shell script used as a wrapper for generating a pairwise Fst table for a given set of populations using vcftools.
# Takes one argument (output file), with two other inputs I manually set here (see below).
# One is a VCF file with all relevant samples included. I used one with one random SNP/locus and 100% SNPs only.
### here "analyses/carunculatus_allPops_Fst-input.vcf"
# The other input is a directory of text files. Each text files has a population name and names of all samples (names as defined in VCF).
### here "vcftools/pops"
#
# I did this the "lazy" way, and didn't make a command line argument for anything but the output file.
# As a matter of principle, I always go with the exact code used to run original analysis as opposed to adding features post-analysis.
# Future implementations would also take the input file and population directory as input. 
# I plan to add a "good" version to https://github.com/ethangyllenhaal/assorted_scripts
#
# Run like: "sh pairwise_Fst_vcftools.sh /path/to/output.txt"

# Just to fill space, really only needs to be a tab 
printf "%s\t" "Pops" > $1
# Prints population file names across the top row
for pop in `ls vcftools/pops`; do 
	printf "%s\t" $pop >> $1 
done
# Adds a newline
echo >> $1 

# Main for loop for making table, nested loop doing every pairwise comparison of populations defined by user.
# Essentially: "For each population, print its name then calculate and print pairwise for each other population."
# Currently relies on internal consistency in order. I don't think this ever ever not the case, but could be remedied by adding a population list input file.
for pop1 in `ls vcftools/pops`; do
	# prints input file name
	printf "%s\t" $pop1 >> $1
	# internal loop for performing pairwise calculation
	for pop2 in `ls vcftools/pops`; do
		# Only calculates pairwise Fst if populations aren't the same.
		if [ $pop1 != $pop2 ]
		then
			# Uses VCFtools to calculate pairwise Fst. However, the output isn't simply a number, and needs to be processed (this is why the output is used).
			vcftools --vcf analyses/carunculatus_allPops_Fst-input.vcf --weir-fst-pop vcftools/pops/$pop1 --weir-fst-pop vcftools/pops/$pop2 --out vcftools/fst_output
			# Locates and outputs the weighted Fst value for the comparison.
			# Note that this doesn't purge this temporary file afterwards.
			grep 'weighted' ./vcftools/fst_output.log | awk '{printf ("%s\t", $7)}' >> $1
		else
			# Prints a "-" as a placeholder for comparisons of the same population.
			# Note that this could be some other value (e.g. nucleotide diversity, Fis) in an analysis.
			printf "%s\t" "-" >> $1
		fi
	done
	# adds a newline after each "pop1" iteration.
	echo >> $1
done

# Removed the .txt after file names. This can/should be changed if you named files differently (e.g. PopName.pop).
sed -i 's/.txt//g' $1
