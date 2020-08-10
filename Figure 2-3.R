#Figures 2-3
library(betareg)
library(mgcv)
library(ggplot2)

#taxa_group here means birds, mammals, and amphibians, file_name is a matrix(CSV format) with foir columns (Latitude, pi_d,hd_d,D_d), pi_d represent nucleotide diversity, hd_d represent haplotype diversity, D_d is the mean average of the geographic distances of each species falling in that latitudinal band  

taxa_group_cytb<-read.csv(file_names,header = T)

taxa_group_cytb$pi_d_tra<-(taxa_group_cytb$pi_d*(length(taxa_group_cytb$pi_d)-1)+0.5)/length(taxa_group_cytb$pi_d)
taxa_group_cytb$hd_d_tra<-(taxa_group_cytb$hd_d*(length(taxa_group_cytb$hd_d)-1)+0.5)/length(taxa_group_cytb$hd_d)
taxa_group_cytb$D_d_tra<-(taxa_group_cytb$D_d*(length(taxa_group_cytb$D_d)-1)+0.5)/length(taxa_group_cytb$D_d)

taxa_group_cytb_all_pi<-betareg(pi_d_tra~I(lat^2)+lat,weights = D_d_tra, data=taxa_group_cytb)
summary(taxa_group_cytb_all_pi)
taxa_group_cytb_all_hd<-betareg(hd_d_tra~I(lat^2)+lat,weights = D_d_tra, data=taxa_group_cytb)
summary(taxa_group_cytb_all_hd)
# you need to illustrated the value of  a, b,c in "method.args = list(start=c(a, b,c))" based on your model
taxa_group_cytb_hd<-ggplot(data=taxa_group_cytb, aes(x=lat, y=hd_d)) +  
  geom_point()  +  
 
  stat_smooth(method = "nls", formula = "y ~ a*x^2+b*x+c",  method.args = list(start=c(a, b,c)), se = FALSE)+labs(x="Latitude", y="Haplotype Diversity",size=16)


taxa_group_cytb_pi<-ggplot(data=taxa_group_cytb, aes(x=lat, y=pi_d)) +  
  geom_point()  +  
 
  stat_smooth(method = "nls", formula = "y ~ a*x^2+b*x+c",  method.args = list(start=c(a, b,c)), se = FALSE)+labs(x="Latitude", y="Nucleotide Diversity",size=16)

#Figures 2-3
library(betareg)
library(mgcv)
library(ggplot2)

#taxa_group here means birds, mammals, and amphibians, file_name is a matrix(CSV format) with foir columns (Latitude, pi_d,hd_d,D_d), pi_d represent nucleotide diversity, hd_d represent haplotype diversity, D_d is the mean average of the geographic distances of each species falling in that latitudinal band  

taxa_group_coi<-read.csv(file_names,header = T)

taxa_group_coi$pi_d_tra<-(taxa_group_coi$pi_d*(length(taxa_group_coi$pi_d)-1)+0.5)/length(taxa_group_coi$pi_d)
taxa_group_coi$hd_d_tra<-(taxa_group_coi$hd_d*(length(taxa_group_coi$hd_d)-1)+0.5)/length(taxa_group_coi$hd_d)
taxa_group_coi$D_d_tra<-(taxa_group_coi$D_d*(length(taxa_group_coi$D_d)-1)+0.5)/length(taxa_group_coi$D_d)

taxa_group_coi_all_pi<-betareg(pi_d_tra~I(lat^2)+lat,weights = D_d_tra, data=taxa_group_coi)
summary(taxa_group_coi_all_pi)
taxa_group_coi_all_hd<-betareg(hd_d_tra~I(lat^2)+lat,weights = D_d_tra, data=taxa_group_coi)
summary(taxa_group_coi_all_hd)
# you need to illustrated the value of  a, b,c in "method.args = list(start=c(a, b,c))" based on your model
taxa_group_coi_hd<-ggplot(data=taxa_group_coi, aes(x=lat, y=hd_d)) +  
  geom_point()  +  
  
  stat_smooth(method = "nls", formula = "y ~ a*x^2+b*x+c",  method.args = list(start=c(a, b,c)), se = FALSE)+labs(x="Latitude", y="Haplotype Diversity",size=16)


taxa_group_coi_pi<-ggplot(data=taxa_group_coi, aes(x=lat, y=pi_d)) +  
  geom_point()  +  
  stat_smooth(method = "nls", formula = "y ~ a*x^2+b*x+c",  method.args = list(start=c(a, b,c)), se = FALSE)+labs(x="Latitude", y="Nucleotide Diversity",size=16)





