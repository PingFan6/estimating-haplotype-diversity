###############################################
# Codes for the paper:
# An approach for estimating haplotype diversity from the sequences with unequal lengths
# Ping Fan, Jon Fjeldså, Xuan Liu, Yafei Dong, Yongbin Chang, Yanhua Qu, Gang Song,Fumin Lei
#A function to simulate the rangdom intraspecific sequences length and calculate the genetic diversity. The algorithm compares all pairwise combinations of sequences and calculates the proportion of loci that differ between the pair. Returns a DataFrame with the identity of sequences, the lengths of overlap region, the number of degenarate bases(including the missing sites),and the computed pairwise divergence values
#'seq_file': A matrix where the rows are aligned genetic sequences, and columns are nucleotide sites, the first columns is the id_number of the sequences. Basepairs must be coded as 1, 2, 3 or 4, or with a 0 signifying that the locus is a degenarate base or missing base. The "id_file" is the assession or id of each sequences, first column is the id_number, the second column is the assession number of each sequences.
#The code in this file by Ping Fan
###############################################


pair_wise_function <- function(seq_file,id_file){
  
  matrix_con <- read.table(seq_file,header = TRUE)
  id_file <- read.table(id_file,header = TRUE)
  matrix_con <- matrix_con_w[,-1]
  data_result<- data.frame(seq1_id=NA,seq2_id=NA,iter=NA,sub_iter=NA,sequence_length=NA,pair_dif=NA,dege_base=NA,num_per_dp=NA)
  num_sequences_size <- length(matrix_con$V1)
  if (num_sequence_size <= 1){
    }
  else 
  {
    id_con <- seq(num_sequence_size)
    iter <- 1
    for (iter in id_con )
    {
      sub_iter <- num_sequence_size
      while(sub_iter - iter > 0)
      {
        seq1 <- matrix_con[iter,]
        seq2 <- matrix_con[sub_iter,]
        length_org <- length(seq1)
       
        s_min <- floor(length_org* 0.5 )
        sequence_length <- floor(runif(1,min = s_min,max = length_org))
        q_max <- length_org - sequence_length + 1
        #
        q_star <- floor(runif(1,min = 1,max = q_max))
        q_fin <- q_star + sequence_length - 1
        pair_dif <- 0
        dege_base <- 0
        for (i in q_star:q_fin)
        {
          if(seq1[i] != seq2[i])
          {
            if(seq1[i] == 0 | seq2[i] ==0){
              dege_base <- dege_base + 1}
            else {
              pair_dif <- pair_dif + 1
            }
          }
        }
        num_per_dp <- pair_dif / (sequence_length - dege_base)
        for (k in 1:length(id_file$iter)) {
          if (iter == id_file$iter[k]){seq1_id = id_file$seqid[k]}
        }
        for (k in 1:length(id_file$iter)) {
          if (sub_iter == id_file$iter[k]){seq2_id = id_file$seqid[k]}
        }
        data1 <- data.frame(seq1_id,seq2_id,iter,sub_iter,sequence_length,pair_dif,dege_base,num_per_dp)
        data_result <- rbind(data_result,data1)
        sub_iter <- sub_iter -1
      }
      iter <- iter +1
      
    }
    data_result <- data_result[-1,]
    write.csv(data_result,seq_random_result,row.names = FALSE)
  }
  
}
