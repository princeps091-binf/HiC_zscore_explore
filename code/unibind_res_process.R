# Examine Unibind output
library(tidyverse)
library(ggridges)
library(svglite)
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
  
  return(tmp_tbl %>% mutate(breast.sample=grepl("breast|mammary",description)) %>% 
           dplyr::select(collection,pValueLog,oddsRatio,support,cellType,breast.sample))
})

breast_TF_tbl<-base::do.call(bind_rows,breast_TF_l)

breast_subset_tbl<-breast_TF_tbl %>% 
  inner_join(.,breast_TF_tbl %>% group_by(collection) %>% 
               summarise(io=any(breast.sample)) %>% 
               filter(io)
  )

breast_TF_tbl %>% 
  group_by(collection) %>% 
  slice_max(pValueLog) %>% 
  ggplot(.,aes(pValueLog))+
  geom_density()


breast_TF_tbl %>% 
  ggplot(.,aes(pValueLog,color=breast.sample))+
  geom_density()



gg_joy<-breast_subset_tbl %>% 
  group_by(collection) %>% 
  filter(n()>5) %>% 
  ggplot(.,aes(pValueLog,y=collection))+
  geom_density_ridges(alpha=0.4)+
  geom_point(aes(color=breast.sample),position="jitter")
gg_joy
ggsave("~/Documents/multires_bhicect/weeklies/weekly49/img/Unibind_breast_joy.svg",gg_joy)

gg_joy<-breast_subset_tbl %>% 
  mutate(set="ALL") %>% 
  bind_rows(.,breast_subset_tbl %>% filter(breast.sample) %>% 
              mutate(set="breast")) %>% 
  ggplot(.,aes(collection,pValueLog,color=set))+
  geom_boxplot()

breast_TF_tbl %>% 
  group_by(collection) %>% 
  slice_max(pValueLog) %>% 
  ggplot(.,aes(pValueLog,support))+
  geom_point()

breast_TF_tbl %>% 
  group_by(collection) %>% 
  slice_max(pValueLog) %>% 
  arrange(desc(pValueLog))
