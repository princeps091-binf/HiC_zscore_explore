library(tidyverse)
library(parallel)
library(furrr)
library(mgcv)

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
