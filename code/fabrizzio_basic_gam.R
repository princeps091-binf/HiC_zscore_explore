library(vroom)
library(mgcv)
library(tidyverse)
options(scipen = 999999999)
#----------------------------------------------
hic_dat_in<-function(dat_file){
  chr_dat<-vroom("PATH/TO/HiC/DATA",delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%
           # ensure HiC score variable is a numeric
           mutate(X3=as.numeric(X3))%>%
           # Filter self-interaction
           filter(!(is.na(X3)))%>%filter(X1!=X2))
}

compute_chr_res_zscore_fn<-function(dat_file,ncluster){
  # Load desired HiC data
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo)%>%
    # Compute genomic ditance between interacting bins
    mutate(d=abs(X1-X2))%>%
    # Perform log transform on distance and HiC score (better behaved input for GAM-computation)
    mutate(lw=log10(X3),ld=log10(d))
  # Compute GAM
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat,cluster = ncluster)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%
    # add predicted HiC values from GAM
    mutate(pred=pred_vec) %>% 
    # Compute zscore using estimation of variance (hic_gam$sig2) and predicted HiC score (pred) from GAM
    mutate(zscore=(lw-pred)/hic_gam$sig2)
  return(chr_dat)
}
