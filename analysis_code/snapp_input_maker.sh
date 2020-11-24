#!/bin/bash

# By: Ethan Gyllenhaal
# Last updated: 21 November 2020
# Wrapper shell script for converting a VCF file to SNAPP input.
# Uses the script below (vcf2nex.pl) to generate initial haploid (0/1) SNAPP input, which I convert here to the combined diploid (0/1/2) SNAPP input.
# Script link: https://github.com/BEAST2-Dev/SNAPP/tree/master/script
# Note that this is a quick and dirty script, and requires lots of values to be changed if you use this for your own data.

# calls perl script to convert vcf to original nexus format
perl vcf2nex.pl < carunc_100_rand_out_UG.vcf > carunc_SNAPP_haploid.nex &&
# outputs the "data" lines of nexus file
sed -n '6,123p' carunc_SNAPP_haploid.nex > snapp_dip.txt &&
# runs python script to make "diploid" individuals for SNAPP
python3 SNAPP_haploid_to_diploid.py &&
# adds nexus header to new nexus file
sed -n '1,5p' carunc_SNAPP_haploid.nex > carunc_SNAPP.nex &&
# halves number of taxa (note that this has to be manually changed
sed -i -e 's/ntax=118/ntax=59/g' carunc_SNAPP.nex &&
# changes symbols
sed -i -e 's/datatype=binary symbols="01"/datatype=standard symbols="012"/g' carunc_SNAPP.nex &&
# adds output of python script to new nexus
cat snapp_out.txt >> carunc_SNAPP.nex &&
# removes '_1" from names
sed -i -e 's/_1 / /g' carunc_SNAPP.nex &&
# adds end of nexus file
echo ';\nEnd;' >> carunc_SNAPP.nex
