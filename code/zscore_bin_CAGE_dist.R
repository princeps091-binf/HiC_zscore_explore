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
