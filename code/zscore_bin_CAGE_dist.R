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
#-------------------------------------------------------------------
## Perform same cage peak count for every HiC bin observed
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"
## Utility functions
hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat <- read_delim(paste0(dat_file,cl_res,"/",chromo,".txt"), 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.nan(X3)))%>%filter(X1!=X2)%>%mutate(d=abs(X1-X2))%>%mutate(lw=log10(X3),ld=log10(d)))
}
compute_bin_cage_overlap_fn<-function(dat_file,cl_res,chromo,full_cage_Grange){
  chr_bin_dat<-hic_dat_in(dat_file,cl_res,chromo) %>% summarise(bins=unique(c(X1,X2))) %>% mutate(chr=chromo,res=cl_res)
  plan(multisession, workers = 3)
  
  chr_bin_dat<-chr_bin_dat %>% mutate(cage.count=unlist(future_pmap(list(chr,bins,res),function(chr,bins,res){
    bin_GRange<-GRanges(seqnames=chr,
                        ranges = IRanges(start=as.numeric(bins),
                                         end=as.numeric(bins)+res_num[res]-1
                        ))
    return(countOverlaps(bin_GRange,full_cage_Grange))
    
    
  })))
  return(chr_bin_dat)
}
compute_chr_res_zscore_fn<-function(dat_file,cl_res,chromo){
  chr_dat<-hic_dat_in(dat_file,cl_res,chromo) 
  hic_gam<-bam(lw~s(ld,bs = "ad"),data = chr_dat)
  pred_vec<-predict(hic_gam,newdata = chr_dat)
  #Compute zscore and predicted HiC magnitude
  chr_dat<-chr_dat%>%mutate(pred=pred_vec,zscore=(chr_dat$lw-pred_vec)/hic_gam$sig2)
  chr_dat<-chr_dat%>%mutate(dist=1/(zscore + abs(min(zscore)-1)))
  return(chr_dat %>% mutate(res=cl_res,chr=chromo))
}  

chromo<-"chr22"
cl_res<-"5kb"

chr_bin_cage_count_tbl<-compute_bin_cage_overlap_fn(dat_file,cl_res,chromo,full_cage_Grange)
chr_zscore_tbl<-compute_chr_res_zscore_fn(dat_file,cl_res,chromo)
#-------------------------------------------------------------------


chr_res_hic_en_tbl <- chr_bin_cage_count_tbl %>% dplyr::slice(1:5) %>% 
      mutate(hic.en=map_dbl(bins,function(x){
      chr_zscore_tbl %>% filter(res==cl_res) %>% filter(X1 == x | X2 == x) %>% summarise(sum(zscore>0)/n()) %>% unlist
    }))

chr_res_hic_en_tbl %>% 
  ggplot(.,aes(cage.count,hic.en))+geom_point()
