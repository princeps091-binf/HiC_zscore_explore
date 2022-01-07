# Examine the zscore properties of CAGE-containing bins
library(tidyverse)
library(GenomicRanges)
library(furrr)
library(mgcv)
options(scipen = 999999999)

res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set
#-------------------------------------------------------------------
# Build the CAGE GRange
cage_coord_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/CAGE/CAGE_coord_tbl.Rda"
out_file<-"./data/CAGE_H1_gene_GRange.Rda"

cage_tbl<-get(base::load(cage_coord_file))
tmp_obj<-names(mget(base::load(cage_coord_file)))
rm(list=tmp_obj)
rm(tmp_obj)

cage_Grange_fn<-function(cage_hmec_a){
  cage_a<-cage_hmec_a%>%filter(!(is.na(start)))
  full_cage_Grange<-   GRanges(seqnames=cage_a$chr,
                               ranges = IRanges(start=cage_a$start,
                                                end=cage_a$end,
                                                names=paste(cage_a$chr,1:nrow(cage_a),sep='_')
                               ))
  mcols(full_cage_Grange)<-tibble(ENSG=cage_a$ENSG)
  return(full_cage_Grange)
}
full_cage_Grange<-cage_Grange_fn(cage_tbl)
#-------------------------------------------------------------------
# Build the hub GRange
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"
res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
chr_set<-unlist(lapply(strsplit(grep('chr',list.files(res_file),value=T),'_'),'[',1))
res_set<- grep('b$',list.files(dat_file),value=T)

#Load the candidate transcription hubs
cl_tbl_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/CAGE/dagger_mres_fdr_01_multi_cagebin_tbl.Rda"

cl_tbl<-get(base::load(cl_tbl_file))
tmp_obj<-names(mget(base::load(cl_tbl_file)))
rm(list=tmp_obj)
rm(tmp_obj)

chr_set<-unique(cl_tbl$chr)
chr_res_l<-vector('list',length(chr_set))
names(chr_res_l)<-chr_set
for (chromo in chr_set){
  print(chromo)
  base::load(paste0(res_file,chromo,"_spec_res.Rda"))
  chr_cl_tbl<-cl_tbl %>% filter(chr==chromo) %>% mutate(bins=chr_spec_res$cl_member[node])
  plan(multisession, workers = 3)
  
  chr_cl_tbl<-chr_cl_tbl %>% mutate(GRange=future_pmap(list(chr,bins,res),function(chr,bins,res){
    return(GRanges(seqnames=chr,
                   ranges = IRanges(start=as.numeric(bins),
                                    end=as.numeric(bins)+res_num[res]-1
                   )))
    
    
  }))
  chr_res_l[[chromo]]<-chr_cl_tbl 
}

chr_res_tbl<-do.call(bind_rows,chr_res_l)
plan(multisession, workers = 3)

chr_res_tbl<-chr_res_tbl %>% mutate(cage_bin_count=future_map(GRange,function(x){
  tmp_vec<-countOverlaps(x,full_cage_Grange)
  return(tibble(bin=start(x),cage.count=tmp_vec))
}))
#-------------------------------------------------------------
## Collect the corresponding hub interaction data with the zscore transformation
## Utility functions
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}
compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo,max.dist,res_num){
  if(max.dist == res_num[cl_res]){
    chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
    chr_dat<-chr_dat %>% filter(d<=median(chr_dat$d))
    
  }  else{
    
    chr_dat<-hic_dat_in(dat_file,cl_res,chromo) %>% filter(abs(X1-X2)<=max.dist)
    
  }
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}  

## Loop through chromosomes and resolution to compute zscore and extrct hub HiC data
chr_hic_dat_l<-vector('list',length(unique(chr_res_tbl$chr)))
names(chr_hic_dat_l)<-unique(chr_res_tbl$chr)

for (chromo in unique(chr_res_tbl$chr)){
  tmp_res_set<-chr_res_combo %>% filter(chr==chromo) %>% distinct(res) %>% unlist
  tmp_res_l<-vector('list',length(tmp_res_set))
  names(tmp_res_l)<-tmp_res_set
  #load the cluster results
  for (cl_res in tmp_res_set){
    message(chromo," : ",cl_res)
    chr_cl_tbl<-chr_res_tbl %>% filter(chr==chromo & res == cl_res) %>% dplyr::select(chr,res,node,bins)
    max.dist<-chr_cl_tbl %>% mutate(max.dist=map_dbl(bins,function(x){
      diff(range(as.numeric(x)))
    })) %>% summarise(max(max.dist)) %>% unlist
    chr_dat<-compute_chr_res_zscore_fn(dat_file,cl_res,chromo,max.dist,res_num)
    cl_bins<-chr_cl_tbl %>% dplyr::select(bins) %>% unnest(cols = c(bins)) %>% distinct() %>% unlist
    cl_dat<-chr_dat %>% filter(X1 %in% as.numeric(cl_bins) & X2 %in% as.numeric(cl_bins))
    rm(chr_dat)
    
    chr_cl_tbl<-chr_cl_tbl %>% mutate(hic.dat=map(bins,function(x){
      return(cl_dat %>% filter(X1 %in% as.numeric(x) & X2 %in% as.numeric(x)))
    }))
    tmp_res_l[[cl_res]]<-chr_cl_tbl
    rm(chr_cl_tbl,cl_bins,max.dist,chr_cl_tbl)
  }  
  chr_hic_dat_l[[chromo]]<-do.call(bind_rows,tmp_res_l)
  
}
chr_hic_dat_tbl<-do.call(bind_rows,chr_hic_dat_l)
save(chr_hic_dat_tbl,file="./data/hub_hic_tbl.Rda")
save(chr_res_tbl,file="./data/hub_cage_tbl.Rda")
#---------------------------------------------------------------------------------------
load("./data/hub_hic_tbl.Rda")
load("./data/hub_cage_tbl.Rda")
plan(multisession, workers = 3)
chr_res_tbl<-chr_res_tbl %>% mutate(cage.bins.vec=future_map(cage_bin_count,function(x){
  return(x %>% filter(cage.count>0) %>% dplyr::select(bin) %>% unlist)
}))
chr_cage_zscore_tbl<-chr_res_tbl %>% dplyr::select(chr,res,node,cage.bins.vec) %>% left_join(.,chr_hic_dat_tbl) %>%
#  dplyr::slice(1234:1245) %>% 
  mutate(cage.hic.vec=future_pmap(list(cage.bins.vec,hic.dat),function(cage.bins.vec,hic.dat){
    hic.dat %>% filter(X1 %in% as.numeric(cage.bins.vec)& X2 %in% as.numeric(cage.bins.vec)) %>% 
      dplyr::select(zscore) %>% unlist
  }))

hub_cage_zscore_tbl<-chr_cage_zscore_tbl %>% dplyr::select(res,cage.hic.vec) %>% unnest(cols=c(cage.hic.vec)) %>% 
  dplyr::rename(zscore=cage.hic.vec) %>% mutate(set='CAGE')
hub_all_zscore_tbl<-chr_hic_dat_tbl %>% mutate(zscore.vec=future_map(hic.dat,function(x){
  return(x$zscore)
})) %>% dplyr::select(res,zscore.vec)%>% unnest(cols=c(zscore.vec)) %>% 
  dplyr::rename(zscore=zscore.vec) %>% mutate(set='ALL')

hub_cage_zscore_tbl %>% bind_rows(.,hub_all_zscore_tbl) %>% 
  mutate(res=fct_relevel(res,names(res_num))) %>% 
  ggplot(.,aes(zscore,color=set))+
  geom_density()+
  facet_wrap(res~.,scales="free")

