# estimating-haplotype-diversity
A new approach for estimating haplotype diversity from mitochondrial DNA

Ping Fan, Jon Fjeldså, Xuan Liu, Yafei Dong, Yongbin Chang, Yanhua Qu, Gang Song,Fumin Lei

# Description
A method to calculate the haplotype diversity from pairwise results. Note that all functions provided in "Calculate pairwise function.R" and "haplotype_diversity_estimate_from_pair_wise_results.R" should be loaded before calculation.

# Data
The alignment DNA sequences that  stored in FASTA format.

# Examples
# 1. Calculate the haplotype divesity 
 Obtain the pairwise result
run_pairwise_function("example_data.fas")
The output file by "run_pairwise_function" is "example_data_result.csv"
Calculate the genetic dviersity 

GD_function("example_data_result.csv")

The out reuslt 
[1] "Haplotype diversity: 0.4 , Nucleotide diversity: 0.000613496932515336"


# 2. Random length analysis
Obtain the simulate pairwise results.Here we set the threshold = 0.5 (50%), repeat number =10

random_pairwise_function("example_data.fas",0.5,10) 

The output file by "random_pairwise_function" is "example_data_1.csv","example_data_2.csv","example_data_3.csv", ... , "example_data_n.csv", where n is the repeat number

Calculate the genetic dviersity 


GD_function("example_data_1.csv")
GD_function("example_data_2.csv")
...
