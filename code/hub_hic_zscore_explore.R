#Process generic aggregate HiC metric for clusters genome-wide
library(tidyverse)
library(mgcv)
library(parallel)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-----------------------------------------
#Utility functions
##input hic data of choice
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

cl_dist_fn<-function(chr_dat,cl_bin_l){
  fn_env<-environment()
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("chr_dat","cl_bin_l"),envir=fn_env)
  cl_dist<-unique(unlist(parLapply(cl,cl_bin_l,function(x){
    tmp_cl_bin<-as.numeric(x)
    return(chr_dat%>%filter(X1 %in% tmp_cl_bin & X2 %in% tmp_cl_bin) %>% mutate(d=abs(X2-X1)) %>% distinct(d) %>% unlist)
  })))
  stopCluster(cl)
  rm(cl)
  return(cl_dist)
}

cl_hic_en_fn<-function(chr_dat,cl_bin_l){
  fn_env<-environment()
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("chr_dat","cl_bin_l"),envir=fn_env)
  cl_dist<-unlist(parLapply(cl,cl_bin_l,function(x){
    tmp_cl_bin<-as.numeric(x)
    return(chr_dat%>%filter(X1 %in% tmp_cl_bin & X2 %in% tmp_cl_bin) %>%   summarise(hic.en=sum(zscore>0)/n()))
  }))
  stopCluster(cl)
  rm(cl)
  return(tibble(cl=names(cl_bin_l),hic.en=cl_dist))
}

cl_hic_dat_fn<-function(chr_dat,cl_bin_l){
  fn_env<-environment()
  cl<-makeCluster(5)
  clusterEvalQ(cl, {
    library(dplyr)
    print("node ready")
  })
  clusterExport(cl,c("chr_dat","cl_bin_l"),envir=fn_env)
  cl_dist<-parLapply(cl,cl_bin_l,function(x){
    tmp_cl_bin<-as.numeric(x)
    return(chr_dat%>%filter(X1 %in% tmp_cl_bin & X2 %in% tmp_cl_bin) )
  })
  stopCluster(cl)
  rm(cl)
  return(tibble(cl=names(cl_bin_l),hic.dat=cl_dist) %>% unnest(cols = c(hic.dat)))
}

#-----------------------------------------
#Load the considered HiC and BHiCect data
dat_file<-"~/Documents/multires_bhicect/data/HMEC/"
res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
chr_set<-unlist(lapply(strsplit(grep('chr',list.files(res_file),value=T),'_'),'[',1))
res_set<- grep('b$',list.files(dat_file),value=T)
#Load the candidate transcription hubs
load('~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/dagger_mres_fdr_01_multi_cagebin_tbl.Rda')

chromo<-"chr15"
tmp_res<-"5kb"
#Subset corresponding hubs
tmp_hub_tbl<-dagger_mres_tbl %>% filter(chr==chromo & res==tmp_res)
#Load the corresponding BHiCect results
load(paste0(res_file,chromo,'_spec_res.Rda'))
cl_bin_l<-chr_spec_res$cl_member[tmp_hub_tbl$node]
chr_dat<-hic_dat_in(dat_file,tmp_res,chromo)
hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
## Build expectation vector for HiC
pred_vec<-predict(hic_gam,newdata = chr_dat)
chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)

# Collect all the interaction distances observed in hubs
cl_dist_set<-cl_dist_fn(chr_dat,cl_bin_l)
chr_dat %>% filter(d %in% cl_dist_set) %>% 
  ggplot(.,aes(zscore))+geom_density()
chr_dat %>% filter(d %in% cl_dist_set) %>% 
  summarise(hic.en=sum(zscore>0)/n())
cl_hic_en_fn(chr_dat,cl_bin_l)%>% ggplot(.,aes(hic.en))+geom_histogram()

cl_kind<-table(unlist(lapply(lapply(strsplit(tmp_hub_tbl$node,split="_"),"[",1:2),function(x)paste(x,collapse = "_"))))
cl_kind_tbl<-tibble(kind=names(cl_kind),n=as.numeric(cl_kind))
rn_cl_set<-unlist(lapply(1:nrow(cl_kind_tbl),function(i){
  sample(grep(cl_kind_tbl$kind[i],names(chr_spec_res$cl_member),value = T),cl_kind_tbl$n[i])
}))
rn_cl_l<-chr_spec_res$cl_member[rn_cl_set]
cl_hic_en_fn(chr_dat,rn_cl_l) %>% ggplot(.,aes(hic.en))+geom_histogram()

cl_hic_dat<-cl_hic_dat_fn(chr_dat,cl_bin_l)
cl_hic_dat %>% 
  ggplot(.,aes(zscore,group=cl))+geom_density()

rn_cl_set<-unlist(lapply(1:nrow(cl_kind_tbl),function(i){
  sample(grep(cl_kind_tbl$kind[i],names(chr_spec_res$cl_member),value = T),cl_kind_tbl$n[i])
}))
rn_cl_l<-chr_spec_res$cl_member[rn_cl_set]
rn_cl_hic_dat<-cl_hic_dat_fn(chr_dat,rn_cl_l)
rn_cl_hic_dat %>% 
  ggplot(.,aes(zscore,group=cl))+geom_density()
