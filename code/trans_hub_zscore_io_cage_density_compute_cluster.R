library(tidyverse)
library(mgcv)
library(furrr)
library(vroom)
library(GenomicRanges)

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
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%
           filter(!(is.na(X3)))%>%filter(X1!=X2)%>%
           mutate(d=abs(X1-X2))%>%
           mutate(lw=log10(X3),ld=log10(d)))
}

compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,res_num,ncluster){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat,cluster = ncluster)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}

compute_bin_cage_overlap_fn<-function(chr_dat,cl_res,chromo,cage_GRange,res_num){
  chr_bin_dat<-chr_dat %>% summarise(bins=unique(c(X1,X2))) %>% mutate(chr=chromo,res=cl_res)
  
  
  bin_GRange<-GRanges(seqnames=chr_bin_dat$chr,
                      ranges = IRanges(start=as.numeric(chr_bin_dat$bins),
                                       end=as.numeric(chr_bin_dat$bins)+res_num[tmp_res]-1
                      ))
  
  bin_count<-countOverlaps(bin_GRange,cage_GRange)
  chr_bin_dat<-chr_bin_dat %>% mutate(cage.count=bin_count)
  
  return(chr_bin_dat)
}

#-----------------------------------------
candidate_hub_file<-"~/data_transfer/candidate_trans_DAGGER_hub/GM12878_union_trans_res_dagger_tbl.Rda"
res_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/spec_res/"
dat_file<-"/storage/mathelierarea/processed/vipin/group/HiC_data/GM12878/"
CAGE_peak_GRange_file<-"~/data_transfer/CAGE_GRange/CAGE_union_GM12878_Grange.Rda"
#-----------------------------------------
#loop through resolution
top_compound_hub_5kb_tbl<-data_tbl_load_fn(candidate_hub_file)
cage_GRange<-data_tbl_load_fn(CAGE_peak_GRange_file)

tmp_res_set<-unique(top_compound_hub_5kb_tbl$res)

for(tmp_res in tmp_res_set){

  # subset hubs to only ones from considered resolution
  hub_chr<-top_compound_hub_5kb_tbl %>% filter(res==tmp_res) %>% distinct(chr) %>% unlist
  # derive available chromosomes at cosidered resolution given HiC data files present
  chr_set<-unlist(lapply(strsplit(list.files(paste0(dat_file,tmp_res)),split="\\."),'[',1))
  # further subset chromosomes to only considerd those for which at hub is present
  chr_set<-chr_set[chr_set %in% hub_chr]
  
  
  # Loop through chromsomes to extract HiC data in vs out of Hub
  io_dat_l<-lapply(chr_set,function(chromo){
    message(chromo,":",tmp_res)
    
    # Load the clusters found in that chromosome
    chr_spec_res<-data_tbl_load_fn(paste0(res_file,chromo,"_spec_res.Rda"))
    
    tmp_hub<-top_compound_hub_5kb_tbl %>% 
      filter(chr==chromo,res==tmp_res) %>% 
      mutate(bins=chr_spec_res$cl_member[node]) %>% 
      mutate(bins=map(bins,function(x) as.numeric(x)))
    
    # Compute z-score for that chromsome at that resolution
    message(chromo," at ",tmp_res,": GAM computation")
    chr_dat<-compute_chr_res_zscore_fn(dat_file,tmp_res,chromo,res_num,10)
    
    # Subset the HiC data to only consider CAGE-interaction
    ## Identify cage-containing bins
    chr_cage_bin_dat<-compute_bin_cage_overlap_fn(chr_dat,tmp_res,chromo,cage_GRange,res_num) %>% filter(cage.count>0)
    ## Subset HiC data accordingly
    chr_cage_hic_dat<-chr_dat %>% filter(X1 %in% chr_cage_bin_dat$bins & X2 %in% chr_cage_bin_dat$bins)
    
    
    hub_chr_dat<-chr_cage_hic_dat %>% filter(X1 %in% unique(unlist(tmp_hub$bins)) & X2 %in% unique(unlist(tmp_hub$bins)))
    
    plan(multisession, workers = 15)
    message(chromo," at ",tmp_res,": Subset hub interaction")
    tmp_hub<-tmp_hub %>% mutate(hub.feature.zscore=future_map(bins,function(x){
      hub_chr_dat %>% filter(X1 %in% x & X2 %in% x) %>% dplyr::select(chr,res,X1,X2,zscore)
    }))
    
    plan(sequential)
    
    tmp_hub_dat<-tmp_hub %>% dplyr::select(hub.feature.zscore) %>% unnest(cols = c(hub.feature.zscore)) %>% mutate(hub.io="hub")
    
    io_hub_dat<-tmp_hub_dat %>% full_join(.,chr_dat %>% dplyr::select(chr,res,X1,X2,zscore))
    
    io_hub_dat<-io_hub_dat %>% mutate(hub.io=ifelse(is.na(hub.io),"out",hub.io),zscore=as.numeric(zscore)) %>% distinct()
    return(io_hub_dat)
    
  })
  
  tmp_res_tbl<-do.call(bind_rows,io_dat_l)
  gg_tmp<-tmp_res_tbl%>%ggplot(.,aes(zscore,color=hub.io))+geom_density()
  ggsave("~/data_transfer/test.png",gg_tmp)
  
  
  
  
  }