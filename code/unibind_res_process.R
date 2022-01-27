# Examine Unibind output
library(tidyverse)
library(ggridges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------
unibind_folder<-"./data/Unibind/UniBind_oneSetBg_HMEC_50kb_hub_vs_HMEC_enh/"

TF_LOLA_res_files<-grep("^col_",list.files(unibind_folder),value=T)

breast_TF_l<-map(TF_LOLA_res_files,function(TF_file){
  message(TF_file)
  tmp_tbl<-read_delim(file=paste0(unibind_folder,TF_file),col_types = list(
    collection = col_character(),
    cellType = col_character()
  ))
  
  return(tmp_tbl %>% filter(grepl("breast",description)) %>% dplyr::select(collection,pValueLog,oddsRatio,support,cellType))
})

breast_TF_tbl<-base::do.call(bind_rows,breast_TF_l)

breast_TF_tbl %>% 
  group_by(collection) %>% 
  slice_max(pValueLog) %>% 
  ggplot(.,aes(pValueLog))+
  geom_density()


breast_TF_tbl %>% 
  ggplot(.,aes(pValueLog,group=collection))+
  geom_density()+ylim(c(0,2))



breast_TF_tbl %>% 
  left_join(.,breast_TF_tbl %>% 
              group_by(collection) %>% 
              summarise(m=mean(pValueLog)) 
  ) %>% 
  mutate(collection=fct_reorder(collection,m)) %>% 
  group_by(collection) %>% 
  filter(n()>3) %>% 
  ggplot(.,aes(pValueLog,y=collection))+
  geom_point(aes(m,collection))+
  geom_density_ridges(alpha=0.4)

breast_TF_tbl %>% 
  group_by(collection) %>% 
  slice_max(pValueLog) %>% 
  ggplot(.,aes(pValueLog,support))+
  geom_point()

breast_TF_tbl %>% 
  group_by(collection) %>% 
  slice_max(pValueLog) %>% 
  arrange(desc(pValueLog))
