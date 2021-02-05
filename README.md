# foulehaioUCE
Scripts and input files used in phylogeographic study of the genus Foulehaio. If you anything here, please cite us (or contact EFG if you're uncertain who to cite):

Mapel, X. M., Gyllenhaal, E. F., Modak, T. H., DeCicco, L. H., Naikatini, A., Utzurrum, R. B., ... & Andersen, M. J. (2021). Inter-and intra-archipelago dynamics of population structure and gene flow in a Polynesian bird. Molecular phylogenetics and evolution, 156, 107034.

This README describes the scripts and input files used for a project assessing the phylogeography of a widespread island generalist (Foulehaio carunculatus) and other related taxa. The data was collected using Illumina sequence of libraries enriched for tetrapod ultraconserved elements. It was then processed using common pipelines for the data (namely phyluce and a variant of seqcap_pop) before input files were generated for population genomic analyses.

______________________________
----- PopGen Input Files -----
______________________________

*_PCA-input.vcf - Input files for various PCAs run. The first word is the species included in that VCF. The second set of words in F. carunculatus files (i.e. before PCA, separted by "-") describes the populations included. Paired with "foulehaio_PCA.R".

carunculatus_allPops_Fst-input.vcf - Input file for calculating pairwise Fst values for all F. carunculatus populations, paired with "pairwise_Fst_vcftools.sh".

carunculatus_ABBA-BABA-input.tsv.gz - Input file for ABBA/BABA tests for F. carunculatus. Paired with "carunculatus_ABBABABA.R".

carunculatus_noMatuku_SNAPP-input.nex - Input file for SNAPP analysis without Matuku individual. This one has populations at a well-divided level, but not totally island specific.

carunculatus_withMatuku_SNAPP-input - Input file for SNAPP analysis without Matuku individual. Note that this one is less population specific.

foulehaio_Fst_pops.zip - Zipped directory of text files used to delimit populations used in generating pairwise Fst values (see pairwise_Fst_vcftools.sh).

____________________________
---- Morphological Data ----
____________________________

Mayr_measures.csv - Input file with all of Mayr's measurements. This does not include any of the Kansas data. Also has the island sizes and latitudes used in analysis.

combined_measures_[unadjusted/adjusted]_[full/noMat].csv - Measurements with weighted averages of Kansas and Mayr measurements. The full/noMat refers to if Matuku is and isn't included respectively (essentially if it is for the RAxML or SNAPP topology). The unadjusted data does not have Kansas data centered by averages where both groups have high sample sizes. The adjusted data does have that centering, and it what is used in our analyses. Also has the island sizes and latitudes used in analysis. All are included for the sake of open data!

foulehaio_measurement_info.xlsx - Info on each measurement dataset. Essentially describes the basics of what went into the wing measurement csvs.

KU_measurements.xlsx - Data used for Kansas calculations. Has a "WingChord" tab with wing chord measurements, a "HandWingIndex" tab with HWI values, and a "Centering" tab with centering calculations.

_________________________
-------- Scripts --------
_________________________

OVERALL: Note that these scripts have only had comments and variable names changed from the original run. This means that they aren't "optimal" in terms of the actual coding. Several have improvements mentioned in the file itself.

pairwise_Fst_vcftools.sh - Script used to calculate pairwise Fst values from a VCF file (carunculatus_allPops_Fst-input.vcf). Full description is in file. Essentially takes an input VCF and a directory with text files defining each population. It then iterates through all pairwise combinations of populations and calculates weighted Fst using VCFtools. The output is a tab-delimited table. Pops used are in "foulehaio_Fst_pops.zip."

carunculatus_PCA.R - Script for generating genomic PCAs to assess population structure in F. carunculatus. This takes in VCF files (the carunculatus_*_PCA-input.vcf series), uses the vcfR package to process input, and finally generates the PCAs using the package adegenet's glPca method. This includes PCAs from several combinations of populations, in order to assess fine scale population structure in the species. Note that this generates PCAs individuals, and manually assigns individuals by VCF order.

taviunensis_PCA.R - Script for generating genomic PCAs to assess population structure in F. taviunensis. This takes in VCF file taviunensis_PCA-input.vcf, uses the vcfR package to process input, and finally generates the PCAs using the package adegenet's glPca method.

procerior_PCA.R - Script for generating genomic PCAs to assess population structure in F. procerior. This takes in VCF file procerior_PCA-input.vcf, uses the vcfR package to process input, and finally generates the PCAs using the package adegenet's glPca method.

carunculatus_ABBABABA.R - Script for running ABBA/BABA analyses on F. carunculatus. This takes in a gzipped, tab delimited format of derived allele freqs made using a script from https://github.com/simonhmartin/genomics_general. Does all combinations in paper, using bootstrap_functions.R to actually run the calculations.

bootstrap_functions.R - Script for running ABBA/BABA bootstraps. This contains f and D statistics, but only D is used here. It generates the value of the statistics, a range of bootstrapped values, bootstrapped p-value, Z value, and a density graph of the values of the bootstraps.

snapp_input_maker.sh - Driver script for converting a VCF to diploid SNAPP input (i.e. Nexus input for processing in Beauti). Converts from VCF to haploid SNAPP using https://github.com/BEAST2-Dev/SNAPP/tree/master/script, then to diploid SNAPP using a custom python script (SNAPP_haploid_to_diploid.py), and finally adding lines for a new nexus file using sed commands. Requires a lot of manual file/value entering by user.

SNAPP_haploid_to_diploid.py - Script for converting a tab-delimited set of haploid nexus SNP data to diploid data. Takes and outputs "middle" lines, header and footer done by "snapp_input_maker.sh".

foulehaio_morphometrics.R - Script used for processing morphometric data in phylogenetic framework.

____________________________
-------- Tree Files --------
____________________________

foulehaio_raxml_tree.tre - Phylogenetic tree of all samples, which was pruned for use in phylogenetic comparative analyses. It is also the phylogeny in Figure 1.

foulehaio_island_snapp_tree.tre - Summary tree of SNAPP posteriors used for phylogenetic comparative analyses, also the "condensed" version of the one shown in Figure 2.

