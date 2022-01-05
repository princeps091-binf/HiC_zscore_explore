library(tidyverse)
library(parallel)
library(furrr)
library(mgcv)
library(Matrix)
library(viridis)
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
cl_res<-"5kb"

chr_dat<-hic_dat_in(dat_file,cl_res,chromo)
hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
pred_vec<-predict(hic_gam,newdata = chr_dat)
#Compute zscore and predicted HiC magnitude
chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)

chr_bin<-unique(c(chr_dat$X1,chr_dat$X2))
options(future.globals.maxSize = 1e9)

plan(multisession, workers = 3)

bin_zscore_tbl<-do.call(bind_rows,future_map(chr_bin,function(bin){
  chr_dat %>% filter(X1==bin | X2==bin) %>% dplyr::select(zscore) %>% mutate(bin=bin)
}))

bin_zscore_tbl %>% group_by(bin) %>% 
  summarise(n=n(),med=min(zscore)) %>% 
  ggplot(.,aes(n,med))+geom_point()

bin_zscore_tbl %>% 
  ggplot(.,aes(zscore,group=bin))+geom_density(line.size=0.1)
