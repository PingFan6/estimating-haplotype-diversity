######################################################################## 
# Codes for the paper: A new approach for estimating haplotype diversity from mitochondrial DNA
#
# Ping Fan, Jon Fjeldså, Xuan Liu, Yafei Dong, Yongbin Chang, Yanhua Qu, Gang Song*, Fumin Lei*
#
# Code in this file by Ping Fan
#########################################################################
# "run_pairwise_function" is a  function to compute the pairwise results from a set of sequences. The algorithm compares all pairwise combinations of sequences and calculates the proportion of loci that differ between the pair. Returns a DataFrame with the identity of sequences, the lengths of each sequence, the overlap, and the computed pairwise divergence values.
# "random_pairwise_function" is a  function to simulate the different overlap range/length of two random chosen sequences with a threshold (e.g. threshold=0.5 means the sequences overlapped by at least 50% of the original length),and compute the pairwise results from the simulation. The algorithm compares all pairwise combinations of sequences and calculates the proportion of loci that differ between the pair. Returns a DataFrame with the identity of sequences, the overlap length, and the computed pairwise divergence values.

#Load the necessary library and functions
library(Biostrings)
library(stringr)

trans_fasta <- function( fasta_name ){
  seq <- readLines(fasta_name) 
  seq <- seq[seq != ""] 
  is.anno <- regexpr("^>", seq, perl=T) 
  seq.anno <- seq[ which(is.anno == 1) ] 
  seq.content <- seq[ which(is.anno == -1) ] 
  start <- which(is.anno == 1) 
  if(length(start) > 1) {
    end <- start[ 2:length(start) ]-1 
    end <- c(end, length(seq) ) 
    distance <- end - start 
    
  index <- 1:length(start) 
  index <- rep(index, distance)
  seqs <- tapply(seq.content, index, paste, collapse="") 
  seq.content<-as.character( seqs ) 
  seq.len <- nchar(seq.content)
  seq.ID <- gsub("^>(\\w+\\|){3}([A-Za-z0-9.]+)\\|.*", "\\2", seq.anno, perl = T) 
  result <- data.frame( seq.ID, seq.anno, seq.len, seq.content ) 
  result
  file_name<-gsub(".fas","_fa.csv",fasta_name)
  write.csv(result,file_name)
  }
}
#Output a matrix where the rows are aligned genetic sequences, and columns are loci. Basepairs are coded as 1, 2, 3 or 4, or with a 0 signifying that the locus is absent or degenerate base from the alignment. Note: the first column is the sequences order marked as: 1,2,3,..., n, where n is the sample size of that species

fasta_matrix <-function(file_name){
    
    library(stringr)
    library(Biostrings)
    test_chose<-read.csv(file_name,header = TRUE)
    test_id_index <-data.frame(iter=test_chose$X,seqid=test_chose$seq.ID)
    test_id_index$seqid1<- gsub(">","",test_id_index$seqid)
   #test_id_index$seqid1<-str_sub(test_id_index$seqid1,1,8)
   test_id_index$seqid<-str_sub(test_id_index$seqid1,1,20)
   judge_value<-"-"
  for (i in 1:length(test_id_index$seqid)) {
    if(grepl(judge_value,test_id_index$seqid[i]))
      {
      
      for (j in 1:str_length(test_id_index$seqid[i])) {
         if(substr(test_id_index$seqid[i],j,j)=="-")
         {
         A_index<-j
       }

      }
      A_index<-A_index + 2
      test_id_index$seqid[i]<-str_sub(test_id_index$seqid[i],1,A_index)  
    }
    else{
      for (j in 1:str_length(test_id_index$seqid[i])) {
      if(substr(test_id_index$seqid[i],j,j)==".")
      {
        A_index<-j
      }
      }
      A_index<-A_index - 1
      test_id_index$seqid[i]<-str_sub(test_id_index$seqid[i],1,A_index)  
    
  } 
  }     
   file_name_index<-gsub("_fa.csv","_index.csv",file_name)

   write.table(test_id_index,file_name_index, row.names = FALSE, na = "",  sep = ",")
                             
    matrix_con <- matrix(ncol = str_length(as.character(test_chose$seq.content[1])),nrow = length(test_chose$seq.ID) )
    for (j in 1:length(test_chose$seq.ID)) {
      
      t_str<-as.character(test_chose$seq.content[j])
      
      #A=1,T=2,G=3,c=4
      for (i in 1:str_length(t_str)) {
        if(substr(t_str,i,i) == "A"){matrix_con[j,i]=1}
        else if(substr(t_str,i,i) == "T"){matrix_con[j,i]=2}
        else if(substr(t_str,i,i) == "G"){matrix_con[j,i]=3}
        else if(substr(t_str,i,i) == "C"){matrix_con[j,i]=4}
        else {matrix_con[j,i]=0}
      }
    }
    matrix_con_w<-as.data.frame(matrix_con)
    file_name_matrix<-gsub("_fa.csv","_matrix.csv",file_name)
    write.csv(matrix_con_w,file_name_matrix)
  }

# "pair_function" use to calculate the pairwise results from matrix file
  
pair_function <-function(file_name){
  matrix_con_w <-read.csv(file_name,header = TRUE)
  id_file_name<- gsub("matrix","index",file_name)
  id_file <-read.csv(id_file_name,header = TRUE)
  matrix_con <-matrix_con_w[,-1]
  data_result<- data.frame(seq1_id=NA,seq2_id=NA,iter=NA,sub_iter=NA,length_seq1=NA,length_seq2=NA,length_common=NA,pair_dif=NA,dege_base=NA,overlap=NA,num_per_dp=NA)
  num_sequence_size <- length(matrix_con$V1)
  if (num_sequence_size <=1){
    data_result_seq1<-data_result
    write.csv(data_result_seq1,file=paste0("seq1","_",file_name,".csv"))}
  else {
    id_con <- seq(num_sequence_size)
    iter <- 1
    for (iter in id_con )
      {
      sub_iter <- num_sequence_size
        while(sub_iter - iter > 0)
        {
          seq1<-matrix_con[iter,]
          seq2<-matrix_con[sub_iter,]
          left_common1 <- min(which(seq1 != 0))
          right_common1 <- max(which(seq1 != 0))
          length_seq1 <-right_common1 - left_common1 +1
          left_common2 <- min(which(seq2 != 0))
          right_common2 <- max(which(seq2 != 0))
          length_seq2 <-right_common2 - left_common2 +1
          left_common <- max(left_common1,left_common2)
          right_common <- min(right_common1,right_common2)
          length_common<- right_common - left_common + 1
          overlap <- length_common/max(length_seq1,length_seq2)
          pair_dif <- 0
          dege_base<-0
          for (i in left_common:right_common)
          {
            if(seq1[i] != seq2[i])
            {
              if(seq1[i] == 0 | seq2[i] ==0){
                dege_base <- dege_base +1}
              else 
                pair_dif <- pair_dif +1
            }
          }
          num_per_dp <- pair_dif / (length_common - dege_base)
          for (k in 1:length(id_file$iter)) {
            if (iter == id_file$iter[k]){seq1_id=id_file$seqid[k]}
          }
          for (k in 1:length(id_file$iter)) {
            if (sub_iter == id_file$iter[k]){seq2_id=id_file$seqid[k]}
          }
          data1<- data.frame(seq1_id,seq2_id,iter,sub_iter,length_seq1,length_seq2,length_common,pair_dif,dege_base,overlap,num_per_dp)
          data_result <- rbind(data_result,data1)
          sub_iter <- sub_iter -1
        }
      progress_status<-"pair_wise_ is running, now finished data %"
      progress<-iter / num_sequence_size * 100
      progress_status<-gsub("data",progress,progress_status)
      print(progress_status)
      iter <- iter +1
    
    }
    data_result<-data_result[-1,]
    file_name_result<-gsub("matrix","result",file_name)
    write.csv(data_result,file_name_result,row.names = FALSE)
    }
  }

# "random_simulate_function" use to simulate the different overlap range/ length with a threshold of 50% and calculate the pairwise results from matrix file 
  
random_simulate_function <-function(file_name,threshold, ith_repeat){
    
    matrix_con_w <-read.csv(file_name,header = TRUE)
    id_file_name<- gsub("matrix","index",file_name)
    id_file <-read.csv(id_file_name,header = TRUE)
    matrix_con <-matrix_con_w[,-1]
    data_result<- data.frame(seq1_id=NA,seq2_id=NA,iter=NA,sub_iter=NA,sequence_length=NA,pair_dif=NA,dege_base=NA,num_per_dp=NA)

    num_sequence_size <- length(matrix_con$V1)
    if (num_sequence_size <=1){
      data_result_seq1<-data_result
      write.csv(data_result_seq1,file=paste0("seq1","_",file_name,".csv"))}
    else 
    {

      id_con <- seq(num_sequence_size)
      iter <- 1
      for (iter in id_con )
      {
        sub_iter <- num_sequence_size
        while(sub_iter - iter > 0)
        {
          seq1<-matrix_con[iter,]
          seq2<-matrix_con[sub_iter,]

          length_org<-length(seq1)
 
          s_min<-floor(length_org*threshold)
          sequence_length<-floor(runif(1,min = s_min,max = length_org))

          q_max<-length_org-sequence_length + 1
          #
          q_star<-floor(runif(1,min = 1,max = q_max))
          q_fin<-q_star+sequence_length-1
          

          pair_dif <- 0

          dege_base<-0
          for (i in q_star:q_fin)
          {
            if(seq1[i] != seq2[i])
            {
              if(seq1[i] == 0 | seq2[i] ==0){
                dege_base <- dege_base +1}
              else {
                pair_dif <- pair_dif +1
              }
            }
          }

          num_per_dp <- pair_dif / (sequence_length - dege_base)
          for (k in 1:length(id_file$iter)) {
            if (iter == id_file$iter[k]){seq1_id=id_file$seqid[k]}
          }
          for (k in 1:length(id_file$iter)) {
            if (sub_iter == id_file$iter[k]){seq2_id=id_file$seqid[k]}
          }
          data1<- data.frame(seq1_id,seq2_id,iter,sub_iter,sequence_length,pair_dif,dege_base,num_per_dp)
          data_result <- rbind(data_result,data1)
          sub_iter <- sub_iter -1
        }
        iter <- iter +1
        
      }
 
      data_result<-data_result[-1,]
      file_name_0<-gsub("_matrix.csv","",file_name)
      
      file_name_result<-paste0(file_name_0,"_",ith_repeat,".csv")
      write.csv(data_result,file_name_result,row.names = FALSE)
    }
    
  }
  
  
  
  
  
  
  
  
# The finally output file of run_pairwise_function is a csv file, whose name contain the "_result" 
run_pairwise_function <- function(fasta_name){
    trans_fasta(fasta_name)
    file.names <- dir(pattern="_fa.csv")
    fasta_matrix(file.names[1])
    file.names <- dir(pattern="_matrix.csv")
    pair_wise_function(file.names[1])
}


# The finally output file of random_pairwise_function are csv file, whose name contain the "_result_i.csv" ,where i is the ith repeat in our simulation.



random_pairwise_fun <- function(fasta_name, threshold, repeat_number){
  trans_fasta(fasta_name)
  file.names <- dir(pattern="_fa.csv")
  fasta_matrix(file.names[1])
  file.names <- dir(pattern="_matrix.csv")
    for (i in 1:repeat_number) {
      
      random_simulate_function(file.names[1],threshold,i)
      progress_status<-"inner_function is running, now finished data %"
      progress<-i 
      progress_status<-gsub("data",progress,progress_status)
      print(progress_status)
      
    }
   
}

