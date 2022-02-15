library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(furrr)
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}



CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_HMEC_Grange.Rda"

cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)

plan(multisession, workers = 3)

enh_H3K27ac_pval_tbl<-do.call(bind_rows,future_map(1:length(cage_enh_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = "~/Documents/multires_bhicect/data/epi_data/HMEC/ENCFF305KUC_pvalue_HMEC_H3K27ac.bigWig",which=cage_enh_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% 
           distinct %>% as_tibble %>% 
           mutate(enh=paste0(seqnames(cage_enh_GRange[x]),"_",start(cage_enh_GRange[x]),"_",end(cage_enh_GRange[x]))))
}))

enh_H3K27ac_fc_tbl<-do.call(bind_rows,future_map(1:length(cage_enh_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = "~/Documents/multires_bhicect/data/epi_data/HMEC/ENCFF981WTU_FC_HMEC_H3K27ac.bigWig",which=cage_enh_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% 
           distinct %>% as_tibble %>% 
           mutate(enh=paste0(seqnames(cage_enh_GRange[x]),"_",start(cage_enh_GRange[x]),"_",end(cage_enh_GRange[x]))))
}))

enh_H3K27ac_tbl<-enh_H3K27ac_fc_tbl %>% dplyr::select(enh,score) %>% full_join(.,enh_H3K27ac_pval_tbl,by=c("enh")) %>% 
  dplyr::rename(fc=score.x,pval=score.y) %>% relocate(enh,seqnames,start,end,strand,width,fc,pval)

  
save(enh_H3K27ac_tbl,file="./data/HMEC_enh_H3K27ac_pval_tbl.Rda")

