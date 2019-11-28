###############################################
# A new approach for estimating haplotype diversity from mitochondrial DNA 
## Code in this file by Ping Fan
###############################################

library(dplyr)
pair_wise_result <- read.data(PW_result,header = TRUE)
sel_result <- filter(pair_wise_result,num_per_dp > 0)
hd <- length(sel_result$num_per_dp)/length(pair_wise_result$num_per_dp)
pi <- mean(pair_wise_result$num_per_dp)
   

 
