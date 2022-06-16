library(tidyverse)
library(mgcv)
library(furrr)
library(vroom)
library(data.tree)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
#Utils. Fn

data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,res_num){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat,cluster = 15)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}  

#-----------------------------------------
candidate_hub_file<-"~/data_transfer/candidate_trans_DAGGER_hub/GM12878_union_trans_res_dagger_tbl.Rda"
res_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
dat_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/"

candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_trans_res_dagger_tbl.Rda"
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"

#-----------------------------------------
options(future.globals.maxSize=20971520000)

trans_hub_tbl<-data_tbl_load_fn(candidate_hub_file)

tmp_chr_set<-unique(trans_hub_tbl$chr)

for(chromo in tmp_chr_set){
  
  chr_spec_res<-data_tbl_load_fn(paste0(res_file,chromo,"_spec_res.Rda"))
  chr_bpt<-FromListSimple(chr_spec_res$part_tree)
  node_ancestor<-chr_bpt$Get(function(x){x$Get('name',traversal='ancestor')})
  node_ancestor<-lapply(node_ancestor,'[',-1)
  
  tmp_hub<-trans_hub_tbl %>%
    filter(chr==chromo) %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>% 
    mutate(bins=map(bins,function(x) as.numeric(x)))
  hub_node<-unlist(tmp_hub%>%dplyr::select(node))
  
  out_hub_cl<-names(node_ancestor)[which(unlist(lapply(node_ancestor,function(x){
    !(any(x %in% hub_node))
  })))]
  #Loop through resolution to extract for eah clusters corresponding HiC data
  chr_res_set<-unique(str_split_fixed(hub_node,"_",2)[,1])
  map(chr_res_set,function(tmp_res){
    hic_dat<-compute_chr_res_zscore_fn(dat_file,tmp_res,chromo,res_num)
    
    tmp_res_hub<-grep(paste0("^",tmp_res,"_"),hub_node,value = T)
    tmp_res_out_cl<-grep(paste0("^",tmp_res,"_"),out_hub_cl,value = T)
    
    tmp_res_hub_dat<-tibble(cl=tmp_res_hub,bins=chr_spec_res$cl_member[tmp_res_hub])
    tmp_res_out_cl_dat<-tibble(cl=tmp_res_out_cl,bins=chr_spec_res$cl_member[tmp_res_out_cl])
    
    tmp_hub_hic_dat<-hic_dat %>% 
      filter(X1 %in% unique(unlist(tmp_res_hub_dat$bins)) & X2 %in% unique(unlist(tmp_res_hub_dat$bins)))
    tmp_out_hic_dat<-hic_dat %>% 
      filter(X1 %in% unique(unlist(tmp_res_out_cl_dat$bins)) & X2 %in% unique(unlist(tmp_res_out_cl_dat$bins)))
    
    plan(multisession,workers=10)
    tmp_res_hub_dat_tbl<-do.call(bind_rows,future_map(tmp_res_hub_dat$bins,function(x){
        return(tmp_hub_hic_dat %>% filter(X1 %in% x & X2 %in% x))
      })) %>% distinct %>% 
      dplyr::select(chr,res,X1,X2,zscore) %>% mutate(io.hub="hub")
    plan(sequential)
    
    plan(multisession,workers=10)
    tmp_res_cl_out_dat_tbl<-do.call(bind_rows,future_map(tmp_res_out_cl_dat$bins,function(x){
      return(tmp_out_hic_dat %>% filter(X1 %in% x & X2 %in% x))
    })) %>% distinct %>% 
      dplyr::select(chr,res,X1,X2,zscore)
    plan(sequential)
    
    tmp_res_dat_tbl<-tmp_res_hub_dat_tbl %>% 
      full_join(.,tmp_res_cl_out_dat_tbl) %>% 
      mutate(io.hub=ifelse(is.na(io.hub),"out",io.hub))
    gg_tmp<-tmp_res_dat_tbl %>% 
      ggplot(.,aes(zscore,color=io.hub))+geom_density()
  })
}