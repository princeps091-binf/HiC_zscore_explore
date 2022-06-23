library(tidyverse)

options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

#-----------------------------------------
res_folder<-"~/data_transfer/candidate_trans_DAGGER_hub/io_hub_out_tbl/GM12878/"
res_file<-"GM12878_io_zscore_"
out_file<-"hub_io_zscore_density_viz_"

tmp_res_set<-str_split_fixed(list.files(res_folder),"_",5)[,4]

for (tmp_res in tmp_res_set){
  message(tmp_res)
  
  tmp_res_tbl<-data_tbl_load_fn(paste0(res_folder,res_file,tmp_res,"_tbl.Rda"))
  
  gg_tmp<-tmp_res_tbl %>% 
    ggplot(.,aes(zscore,color=hub.io))+
    geom_density()+
    scale_color_brewer(palette="Set1")+
    theme_minimal()+
    facet_wrap(res~.)
  
  ggsave(paste0(res_folder,out_file,tmp_res,".svg"),gg_tmp)
  
}
