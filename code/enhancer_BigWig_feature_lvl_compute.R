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



CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_GM12878_Grange.Rda"
pval_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/CTCF/ENCODE/BigWig/ENCFF749HDD_GM12878_PVAL.bigWig"
fc_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/CTCF/ENCODE/BigWig/ENCFF312KXX_GM12878_FC.bigWig"

out_file<-"./data/GM12878_enh_CTCF_lvl_tbl.Rda"

cage_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)

plan(multisession, workers = 3)

enh_feature_pval_tbl<-do.call(bind_rows,future_map(1:length(cage_enh_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = pval_file,which=cage_enh_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% 
           distinct %>% as_tibble %>% 
           mutate(enh=paste0(seqnames(cage_enh_GRange[x]),"_",start(cage_enh_GRange[x]),"_",end(cage_enh_GRange[x]))))
}))

enh_feature_fc_tbl<-do.call(bind_rows,future_map(1:length(cage_enh_GRange),function(x){
  H3K27ac_hub_enh<-rtracklayer::import(con = fc_file,which=cage_enh_GRange[x])
  return(as.data.frame(H3K27ac_hub_enh) %>% 
           distinct %>% as_tibble %>% 
           mutate(enh=paste0(seqnames(cage_enh_GRange[x]),"_",start(cage_enh_GRange[x]),"_",end(cage_enh_GRange[x]))))
}))

enh_feature_tbl<-enh_feature_fc_tbl %>% dplyr::select(enh,score) %>% full_join(.,enh_feature_pval_tbl,by=c("enh")) %>% 
  dplyr::rename(fc=score.x,pval=score.y) %>% relocate(enh,seqnames,start,end,strand,width,fc,pval)


save(enh_feature_tbl,file=out_file)

