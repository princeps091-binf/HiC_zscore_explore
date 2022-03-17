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



CAGE_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_tss_HMEC_Grange.Rda"
pval_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/ENCFF305KUC_pvalue_HMEC_H3K27ac.bigWig"
fc_file<-"~/Documents/multires_bhicect/data/epi_data/HMEC/ENCFF981WTU_FC_HMEC_H3K27ac.bigWig"

out_file<-"./data/HMEC_tss_H3K27ac_lvl_tbl.Rda"

cage_GRange<-data_tbl_load_fn(CAGE_GRange_file)

plan(multisession, workers = 5)

enh_feature_pval_tbl<-do.call(bind_rows,future_map(1:length(cage_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = pval_file,which=cage_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% 
           distinct %>% as_tibble %>% 
           mutate(enh=paste0(seqnames(cage_GRange[x]),"_",start(cage_GRange[x]),"_",end(cage_GRange[x]))))
}))

enh_feature_fc_tbl<-do.call(bind_rows,future_map(1:length(cage_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = fc_file,which=cage_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% 
           distinct %>% as_tibble %>% 
           mutate(enh=paste0(seqnames(cage_GRange[x]),"_",start(cage_GRange[x]),"_",end(cage_GRange[x]))))
}))

enh_feature_tbl<-enh_feature_fc_tbl %>% dplyr::select(enh,score) %>% full_join(.,enh_feature_pval_tbl,by=c("enh")) %>% 
  dplyr::rename(fc=score.x,pval=score.y) %>% relocate(enh,seqnames,start,end,strand,width,fc,pval)


save(enh_feature_tbl,file=out_file)

