library(tidyverse)
library(GenomicRanges)
library(furrr)
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
Build_coord_fn<-function(top_compound_hub_5kb_tbl,spec_res_file){
  coord_tbl<-do.call(bind_rows,map(unique(top_compound_hub_5kb_tbl$chr),function(chromo){
    message(chromo)
    base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
    tmp_tbl<-top_compound_hub_5kb_tbl %>% 
      filter(chr==chromo) %>% 
      mutate(bins=chr_spec_res$cl_member[parent.hub]) %>% 
      mutate(bins=map(bins,as.numeric)) 
    
  }))
}

Build_GRange_fn<-function(chromo,res,bins,res_num){
  inter_cl_Grange<-   GRanges(seqnames=chromo,
                              ranges = IRanges(start=bins,
                                               end=bins + res_num[res]-1
                              ))
  inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
  return(inter_cl_Grange)
  
}

#-----------------------------------------
candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/candidate_compound_hub/HMEC_5kb_tss_compound_hub.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/HMEC/spec_res/"

CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_HMEC_Grange.Rda"
CAGE_enh_H3K27ac_file<-"./data/HMEC_enh_H3K27ac_pval_tbl.Rda"
#-----------------------------------------

compound_hub_5kb_tbl<-data_tbl_load_fn(candidate_hub_file)

top_compound_hub_5kb_tbl<-do.call(bind_rows,map(unique(compound_hub_5kb_tbl$chr),function(chromo){
  tmp_tbl<-compound_hub_5kb_tbl %>% 
    filter(chr==chromo)
  return(tmp_tbl %>% 
           filter(parent.hub %in% tmp_tbl$parent.hub[which(!(tmp_tbl$parent.hub %in% unique(unlist(tmp_tbl$ch.hub))))])
  )
  
}))

top_compound_hub_5kb_tbl<-Build_coord_fn(top_compound_hub_5kb_tbl,spec_res_file) %>% 
  mutate(GRange=pmap(list(chr,res,bins),function(chromo,res,bins){
    Build_GRange_fn(chromo,res,bins,res_num)
  }))


cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)
enh_H3K27ac_tbl<-data_tbl_load_fn(CAGE_enh_H3K27ac_file)


top_compound_hub_5kb_tbl<-top_compound_hub_5kb_tbl %>% 
  mutate(peak.content=map(GRange,function(x){
    tmp_hit<-unique(queryHits(findOverlaps(cage_enh_GRange,x)))
    paste(seqnames(cage_enh_GRange)[tmp_hit],start(cage_enh_GRange)[tmp_hit],end(cage_enh_GRange)[tmp_hit],sep="_")
  }))

in_set<-top_compound_hub_5kb_tbl %>% 
#  filter(res=="5kb") %>% 
  dplyr::select(peak.content) %>% 
  unnest(cols=c(peak.content)) %>% 
  distinct() %>% unlist

enh_H3K27ac_tbl %>% 
  mutate(hub.io=ifelse(enh %in% in_set,"in","out")) %>% 
  ggplot(.,aes(pval,color=hub.io))+geom_density()+scale_x_log10()+xlab("-log(p-value)")

zero_tresh<-10**(floor(min(log10(enh_H3K27ac_tbl$fc[enh_H3K27ac_tbl$fc>0]))) -1)
enh_H3K27ac_tbl %>% 
  mutate(hub.io=ifelse(enh %in% in_set,"in","out")) %>% 
  ggplot(.,aes(fc+zero_tresh,color=hub.io))+geom_density()+
  scale_x_log10(breaks=c(zero_tresh,0.1,1,10,100),labels=c(0,0.1,1,10,100))+ 
  xlab("fold-change")

