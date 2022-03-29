library(GenomicRanges)
library(tidyverse)

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



tbl_files<-grep("*_io_tbl.Rda",list.files("./data/hub_io_zscore/top_hub/CTCF/HMEC/"),value = T)
io_hub_dat_tbl<-lapply(paste0("./data/hub_io_zscore/top_hub/CTCF/HMEC/",tbl_files),function(file){
  return(data_tbl_load_fn(file))
})

io_hub_dat_tbl<-do.call(bind_rows,io_hub_dat_tbl)

io_hub_dat_tbl %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(zscore,color=hub.io))+
  geom_density()+facet_wrap(res~.,scales="free")

ggsave("~/Documents/multires_bhicect/weeklies/weekly54/img/CTCF_HMEC_zscore_io_dens.png")
