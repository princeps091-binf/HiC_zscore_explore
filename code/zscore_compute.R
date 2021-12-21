library(tidyverse)
library(parallel)
library(furrr)
library(mgcv)
library(Matrix)
library(viridis)
library(seriation)
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
chr_dat<-chr_dat%>%mutate(dist=1/(zscore + abs(min(zscore)-1)))

chr_dat %>% ggplot(.,aes(dist))+geom_histogram()
chr_dat %>% ggplot(.,aes(dist,zscore))+geom_point()
#-----------------------------------------
#-----------------------------------------
# visualisation of matrix
full_f_mat<-function(cl_mat,res,x_col){
  
  range_5kb<-range(unique(c(cl_mat$X1,cl_mat$X2)))
  bin_5kb<-seq(range_5kb[1],range_5kb[2],by=res)
  #add the bins not present in original Hi-C dataset
  #miss_bin<-bin_5kb[which(!(bin_5kb %in% unique(c(mat_df$X1,mat_df$X2))))]
  
  id_conv<-seq_along(bin_5kb)
  names(id_conv)<-bin_5kb
  
  cl_mat$ego_id<-id_conv[as.character(cl_mat$X1)]
  cl_mat$alter_id<-id_conv[as.character(cl_mat$X2)]
  
  chr_mat<-sparseMatrix(i=cl_mat$ego_id,cl_mat$alter_id,x=as.numeric(cl_mat[[x_col]]),symmetric = T)
  dimnames(chr_mat)<-list(names(id_conv),names(id_conv))
  return(chr_mat)
}

chr_mat<-full_f_mat(chr_dat,res_num[cl_res],"zscore")
image(as.matrix(chr_mat),col=viridis(100))
image(cor(as.matrix(chr_mat)),col=viridis(100))

chr_dist_mat<-as.dist(as.matrix(full_f_mat(chr_dat,res_num[cl_res],"dist")))
chr_dist_mat<-(1-cor(as.matrix(chr_mat)))
chr_dist_mat[is.na(chr_dist_mat)]<-3
chr_dist_mat<-as.dist(chr_dist_mat)
o <- seriate(chr_dist_mat,method = "MDS_metric")
image(as.matrix(chr_mat)[get_order(o),get_order(o)],col=viridis(100))
#-------------------------------
#Produce the eigenvector and corresponding re-ordering 
# => Clear correspondence up til sign inversion with Spectral/MDS-metric re-ordering of bins through seriation package
chr_cor_mat<-cor(as.matrix(chr_mat),use = "na.or.complete")
chr_cor_mat[is.na(chr_cor_mat)]<-0
cor_mat_svd<-svd(chr_cor_mat,nv = 2,nu=2)
eigen_tbl<-tibble(bin=as.numeric(rownames(chr_cor_mat)),eigen1=cor_mat_svd$u[,1]) %>%   
  mutate(ID=1:n()) %>% 
  arrange(eigen1) 
image(as.matrix(chr_mat)[eigen_tbl$ID,eigen_tbl$ID],col=viridis(100))
