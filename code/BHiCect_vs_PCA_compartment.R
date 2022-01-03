library(tidyverse)
library(mgcv)
library(Matrix)
library(viridis)
library(seriation)
library(caret)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
## Examine how BHiCect divides the chromosome

### Compute the z-score
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
rm(pred_vec,hic_gam)
#-----------------------------------------
## Visualisation
### Build the matrix
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
### Build the correlation matrix
chr_cor_mat<-cor(as.matrix(chr_mat),use = "na.or.complete")
#### replace NA with 0-correlations
chr_cor_mat[is.na(chr_cor_mat)]<-0


lp_fn<-function(x){
  
  Dinv=Diagonal(nrow(x),1/Matrix::rowSums(x))
  
  lp_chr1=Diagonal(nrow(x),1)-Dinv %*% x
  
  temp<-eigen(lp_chr1)
  return(list(vectors=as_tibble(temp[['vectors']][,c(length(temp$values)-1,length(temp$values))],.name_repair = "universal"),values=temp[['values']][c(length(temp$values)-1,length(temp$values))]))
    
  }

Spectral_cl_tbl<-lp_fn(1e6+chr_cor_mat)$vectors
Spectral_cl_tbl<-Spectral_cl_tbl %>% 
  rename(eigen2="...1") %>% 
  mutate(ID=1:n()) %>% 
  arrange(eigen2) 

image(chr_cor_mat[Spectral_cl_tbl$ID,Spectral_cl_tbl$ID],col=viridis(100))

cor_mat_svd<-svd(chr_cor_mat,nv = 2,nu=2)
eigen_tbl<-tibble(bin=as.numeric(rownames(chr_cor_mat)),eigen1=cor_mat_svd$u[,1]) %>%   
  mutate(ID=1:n()) %>% 
  arrange(eigen1) 
image(chr_cor_mat[eigen_tbl$ID,eigen_tbl$ID],col=viridis(100))

eigen_tbl %>% inner_join(.,Spectral_cl_tbl) %>% ggplot(.,aes(eigen1,eigen2))+geom_point()
#-----------------------------------------
## Compare with the original spcetral separation done by current BHiCect method
power_trans_fn<-function(x){
  preprocessParams <- BoxCoxTrans(x$X3,na.rm = T)
  x <- x %>% mutate(weight=predict(preprocessParams, x$X3))
  x <- x %>% mutate(weight=weight+(1-min(weight,na.rm = T)))
  return(x)
}
chr_dat<-chr_dat %>% power_trans_fn()
chr_pow_mat<-full_f_mat(chr_dat,res_num[cl_res],"weight")
chr_raw_mat<-full_f_mat(chr_dat,res_num[cl_res],"X3")

BHiCect_cl_tbl<-lp_fn(chr_pow_mat)$vectors
BHiCect_cl_tbl<-BHiCect_cl_tbl %>% 
  rename(eigen2="...1") %>% 
  mutate(ID=1:n()) %>% 
  arrange(eigen2) 

image(as.matrix(chr_pow_mat)[BHiCect_cl_tbl$ID,BHiCect_cl_tbl$ID],col=viridis(100))
