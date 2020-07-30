###############################################
# Codes for the paper:
# A new approach for estimating haplotype diversity from mitochondrial DNA 
# Ping Fan, Jon Fjeldså, Xuan Liu, Yafei Dong, Yongbin Chang, Yanhua Qu, Gang Song,Fumin Lei
## Code in this file by Ping Fan
###############################################
# PW_result is the result output by run_pairwise_function/random_pairwise_function, whose file name contain the "result".
# GD_function is use to calaulate the haplotype diversity and nucleotide diversity from pairwise results. 
library(dplyr)
GD_function <- function(PW_result){
pair_wise_result <- read.csv(PW_result,header = TRUE)
sel_result <- filter(pair_wise_result,num_per_dp > 0)
hd <- length(sel_result$num_per_dp)/length(pair_wise_result$num_per_dp)
pi <- mean(pair_wise_result$num_per_dp)
result_output<-"Haplotype diversity: hd , Nucleotide diversity: pi"
    
      result_output<-gsub("hd",hd,result_output)
      result_output<-gsub("hd",pi,result_output)
      print(result_output)
}
   

 
