library(tidyverse)
library(GenomicRanges)
library(mgcv)
library(furrr)
library(vroom)
library(svglite)
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
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F,show_col_types = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}

compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,res_num){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}  

compute_bin_cage_overlap_fn<-function(chr_dat,cl_res,chromo,full_cage_Grange,res_num){
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
candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/HMEC_union_dagger_tbl.Rda"
CAGE_peak_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_union_HMEC_Grange.Rda"
CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_HMEC_Grange.Rda"
CAGE_enh_tbl_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/CAGE/CAGE_enh_tbl.Rda"


res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"
dat_file<-"~/Documents/multires_bhicect/data/HMEC/"


dagger_hub_tbl<-data_tbl_load_fn(candidate_hub_file)

cage_GRange<-data_tbl_load_fn(CAGE_peak_GRange_file)

cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)

cage_enh_tbl<-data_tbl_load_fn(CAGE_enh_tbl_file)

tmp_res<-"50kb"
hub_chr<-dagger_hub_tbl %>% filter(res==tmp_res) %>% distinct(chr) %>% unlist
chr_set<-unlist(lapply(strsplit(list.files(paste0(dat_file,tmp_res)),split="\\."),'[',1))
chr_set<-chr_set[chr_set %in% hub_chr]

enh_Grange_l<-lapply(chr_set,function(chromo){
  message(chromo)
  chr_spec_res<-data_tbl_load_fn(paste0(res_file,chromo,"_spec_res.Rda"))
  
  chr_dat<-compute_chr_res_zscore_fn(dat_file,tmp_res,chromo,res_num)
  chr_cage_bin_dat<-compute_bin_cage_overlap_fn(chr_dat,tmp_res,chromo,cage_GRange,res_num) %>% filter(cage.count>0)
  chr_cage_hic_dat<-chr_dat %>% filter(X1 %in% chr_cage_bin_dat$bins & X2 %in% chr_cage_bin_dat$bins)
  
  tmp_hub<-dagger_hub_tbl %>% 
    filter(chr==chromo,res==tmp_res) %>% 
    mutate(bins=chr_spec_res$cl_member[node]) %>% 
    mutate(bins=map(bins,function(x) as.numeric(x)))
  
  tmp_hub<-tmp_hub %>% mutate(hub.cage.zscore=future_map(bins,function(x){
    chr_cage_hic_dat %>% filter(X1 %in% x & X2 %in% x) %>% dplyr::select(chr,res,X1,X2,zscore)
  }))
  
  tmp_hub_dat<-tmp_hub %>% dplyr::select(hub.cage.zscore) %>% unnest(cols = c(hub.cage.zscore)) %>% mutate(hub.io="hub")
  
  hub_cage_bin_inter_tbl<-tibble(chr=chromo,res=tmp_res,bin=unique(c(tmp_hub_dat$X1,tmp_hub_dat$X2))) %>% 
    mutate(end= bin +res_num[res]-1)
  
  hub_cage_bin_GRange<-reduce(GRanges(seqnames=hub_cage_bin_inter_tbl$chr,
                                      ranges = IRanges(start=hub_cage_bin_inter_tbl$bin,
                                                       end=hub_cage_bin_inter_tbl$end
                                      )))
  
  return(cage_enh_GRange[unique(queryHits(findOverlaps(cage_enh_GRange,hub_cage_bin_GRange)))])
})
hub_50kb_enh_Grange<-unlist(GRangesList(enh_Grange_l))
plan(multisession, workers = 3)

hmec_enh_H3K27ac_tbl<-do.call(bind_rows,future_map(1:length(cage_enh_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = "~/Documents/multires_bhicect/data/epi_data/HMEC/ENCFF305KUC_pvalue_HMEC_H3K27ac.bigWig",which=cage_enh_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% distinct %>% as_tibble %>% mutate(enh=paste0(seqnames(cage_enh_GRange[x]),"_",start(cage_enh_GRange[x]),"_",end(cage_enh_GRange[x]))))
}))

hmec_enh_H3K27ac_fc_tbl<-do.call(bind_rows,future_map(1:length(cage_enh_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = "~/Documents/multires_bhicect/data/epi_data/HMEC/ENCFF981WTU_FC_HMEC_H3K27ac.bigWig",which=cage_enh_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% distinct %>% as_tibble %>% mutate(enh=paste0(seqnames(cage_enh_GRange[x]),"_",start(cage_enh_GRange[x]),"_",end(cage_enh_GRange[x]))))
}))

hub_enh_tbl<-tibble(hub.enh=paste0(seqnames(hub_50kb_enh_Grange),"_",start(hub_50kb_enh_Grange),"_",end(hub_50kb_enh_Grange)))

hmec_enh_H3K27ac_fc_tbl %>% full_join(.,hmec_enh_H3K27ac_tbl,by=c("enh")) %>% 
  mutate(hub.io=ifelse(enh %in% hub_enh_tbl$hub.enh,"in","out")) %>% 
  ggplot(.,aes(score.y,color=hub.io))+geom_density()+scale_x_log10()+xlab("-log(p-value)")
ggsave("~/Documents/multires_bhicect/weeklies/weekly49/img/hub_io_H3K27ac_pval_enh.svg")
hmec_enh_H3K27ac_fc_tbl %>% full_join(.,hmec_enh_H3K27ac_tbl,by=c("enh")) %>% 
  mutate(hub.io=ifelse(enh %in% hub_enh_tbl$hub.enh,"in","out")) %>% 
  ggplot(.,aes(score.x,color=hub.io))+geom_density()+scale_x_log10()+xlab("fold-change")
ggsave("~/Documents/multires_bhicect/weeklies/weekly49/img/hub_io_H3K27ac_fc_enh.svg")

hub_enh_tbl<-hub_enh_tbl %>% mutate(chr=as.character(seqnames(hub_50kb_enh_Grange)),start=start(hub_50kb_enh_Grange),end=end(hub_50kb_enh_Grange))

cage_enh_tbl %>% left_join(.,hub_enh_tbl) %>% mutate(hub.enh=ifelse(is.na(hub.enh),"out","in")) %>% 
  ggplot(.,aes(m,color=hub.enh))+geom_density()+scale_x_log10()+xlab("CAGE I.")
ggsave("~/Documents/multires_bhicect/weeklies/weekly49/img/hub_io_CAGE_I_enh.svg")

hmec_enh_H3K27ac_tbl %>% dplyr::select(enh,score) %>% 
  full_join(.,cage_enh_tbl %>% mutate(enh=paste(chr,start,end,sep="_"))) %>% 
  mutate(hub.io=ifelse(enh %in% hub_enh_tbl$hub.enh,"in","out")) %>% 
  ggplot(.,aes(m,score))+
  geom_point(alpha=0.2) +
  scale_x_log10()+scale_y_log10()+
  xlab("mean CAGE I")+ylab("H3K27ac -log(pvalue)")+facet_wrap(hub.io~.)
ggsave("~/Documents/multires_bhicect/weeklies/weekly49/img/hub_io_CAGE_I_vs_H3K27ac_pvalue_enh.svg")
