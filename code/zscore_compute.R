library(tidyverse)
library(parallel)
library(furrr)
library(mgcv)
library(Matrix)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
## Produce GAM for zscore for every chromosome at every resolution
##input hic data of choice
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}
dat_file<-"~/Documents/multires_bhicect/data/HMEC/"
chromo<-"chr22"
cl_res<-"100kb"

chr_dat<-hic_dat_in(dat_file,cl_res,chromo)
hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
pred_vec<-predict(hic_gam,newdata = chr_dat)
#Compute zscore and predicted HiC magnitude
chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
#-----------------------------------------
# visualisation of matrix
full_f_mat<-function(cl_mat,res){
  
  range_5kb<-range(unique(c(cl_mat$X1,cl_mat$X2)))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat$zscore),symmetric = T)
  
}

chr_mat<-full_f_mat(chr_dat,res_num[cl_res])
image(as.matrix(chr_mat))
