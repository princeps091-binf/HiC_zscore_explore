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

#-----------------------------------------
TAD_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/pval_tbl/CAGE_union_HMEC_TAD_pval_tbl.Rda"

CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_HMEC_Grange.Rda"
CAGE_enh_H3K27ac_file<-"./data/HMEC_enh_H3K27ac_pval_tbl.Rda"
#-----------------------------------------

TAD_tbl<-data_tbl_load_fn(TAD_file)

CAGE_enh_H3K27ac_tbl<-data_tbl_load_fn(CAGE_enh_H3K27ac_file)
CAGE_enh_GRange<-data_tbl_load_fn(CAGE_enh_GRange_file)
mcols(CAGE_enh_GRange)<-tibble(enh=paste0(seqnames(CAGE_enh_GRange),"_",start(CAGE_enh_GRange),"_",end(CAGE_enh_GRange)))

TAD_set<-TAD_tbl %>%
  group_by(chr) %>% 
  mutate(FDR=p.adjust(emp.pval,method='fdr')) %>% 
  filter(FDR<= 0.01)

TAD_set_GRange<-reduce(do.call("c",TAD_set$GRange))

TAD_in_enh<-mcols(CAGE_enh_GRange)$enh[unique(queryHits(findOverlaps(CAGE_enh_GRange,TAD_set_GRange)))]

zero_tresh<-10**(floor(min(log10(CAGE_enh_H3K27ac_tbl$fc[CAGE_enh_H3K27ac_tbl$fc>0]))) -1)

CAGE_enh_H3K27ac_tbl %>% 
  mutate(tad.io=ifelse(enh %in% TAD_in_enh,"in","out")) %>% 
  ggplot(.,aes(fc+zero_tresh,color=tad.io))+geom_density()+
  scale_x_log10(breaks=c(zero_tresh,0.1,1,10,100),labels=c(0,0.1,1,10,100))+ 
  scale_color_brewer(palette="Set1")+
  xlab("fold-change")
