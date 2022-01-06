library(tidyverse)
library(mgcv)

#-----------------------------------------
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}
#-----------------------------------------
rao_loop <- read_delim("~/Documents/multires_bhicect/data/GM12878/rao_loop.txt",  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#-----------------------------------------
#Load the considered HiC and BHiCect data
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"

chromo<-"chr1"
rao_chromo<-1
tmp_res<-"5kb"

chr_dat<-hic_dat_in(dat_file,tmp_res,chromo)
hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
## Build expectation vector for HiC
pred_vec<-predict(hic_gam,newdata = chr_dat)
chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)

chr_loop_tbl<-rao_loop %>% filter(chr1==rao_chromo) %>% 
  filter(abs(x1-x2)==res_num[tmp_res]) %>% 
  dplyr::select(chr1,x1,y1) %>% 
  mutate(chr=paste0('chr',chr1)) %>% 
  rename(X1=x1,X2=y1) %>% 
  dplyr::select(chr,X1,X2) %>% 
  mutate(loop=T)



chr_dat %>% left_join(.,chr_loop_tbl) %>% 
  filter(loop) %>% 
  ggplot(.,aes(ld,zscore,color=loop))+geom_point()+geom_point(aes(ld,pred),color="black")
