library(tidyverse)
library("MCMCglmm")
library(data.table)
library(pedigree)

setwd("/Volumes/Seagate Portable Drive")
birth_wt_df=read.table("PhD/Chapter_4_inbr_dep/birth_wt_df.txt", sep=",", header=T)

birth_wt_df_na_rm=birth_wt_df%>%
  na.omit()


FROH<-read.table("PhD/2023_ROH_search/2021_sleuthed_052023.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)%>%filter(nchar(Code)==5)

MB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) %>%#completes for all chromosomes 
  mutate(MB_chr=KB_chr/1000)%>%
  select(-KB_chr)%>%
  reshape2::dcast(Code~CHR)%>%
  mutate(ROHsum = rowSums(.[2:34]))%>%mutate(ROH_sum_div=ROHsum/33)

colnames(MB_perLG) <- c("Code", paste0("ROH_chr", 1:33),"ROHsum","ROH_sum_div")

birth_wt_df_na_rm=birth_wt_df_na_rm%>%
  left_join(MB_perLG)


ped=read.table("PhD/Chapter_4_inbr_dep/pedigree.txt", sep=",", header=T)

##sort out pedigree
ids_in<-as.matrix(birth_wt_df_na_rm%>%select(Code))
pruned<-prunePed(ped, ids_in)
ord<-orderPed(pruned)
ped_ordered=pruned[order(ord),]
ped_ordered$Code <- as.factor(ped_ordered$Code)
ped_ordered$MumCode <- as.factor(ped_ordered$MumCode)
ped_ordered$Sire <- as.factor(ped_ordered$Sire)

Ainv<-inverseA(ped_ordered, nodes="ALL")$Ainv

#####################
### making variables factors 

birth_wt_df_na_rm$MotherStatus[birth_wt_df_na_rm$MotherStatus=='Na\xefve']<-'Naive'

birth_wt_df_na_rm$Code <- as.factor(birth_wt_df_na_rm$Code)
birth_wt_df_na_rm$MumCode <- as.factor(birth_wt_df_na_rm$MumCode)
birth_wt_df_na_rm$MotherStatus <- as.factor(birth_wt_df_na_rm$MotherStatus)
birth_wt_df_na_rm$BirthYear <- as.factor(birth_wt_df_na_rm$BirthYear)
birth_wt_df_na_rm$Sex <- as.factor(birth_wt_df_na_rm$Sex)


table(birth_wt_df_na_rm$MotherStatus)

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#### RUNNING MULTI-MEMBERSHIP MODEL ####
birth_wt_model_ROHONLY<-MCMCglmm(CaptureWt~1 + Sex + AgeHrs + MotherStatus + mum_age+mum_age_sq+Day_seq+ROH_sum_div, 
                random= ~ Code +
                  idv(ROH_chr1+ROH_chr2+ROH_chr3+ROH_chr4+ROH_chr5+ROH_chr6+ROH_chr7+ROH_chr8+
                        ROH_chr9+ROH_chr10+ROH_chr11+ROH_chr12+ROH_chr13+ROH_chr14+ROH_chr15+ROH_chr16+
                        ROH_chr17+ROH_chr18+ROH_chr19+ROH_chr20+ROH_chr21+ROH_chr22+ROH_chr23+ROH_chr24+
                        ROH_chr25+ROH_chr26+ROH_chr27+ROH_chr28+ROH_chr29+ROH_chr30+ROH_chr31+ROH_chr32+
                        ROH_chr33)+BirthYear +MumCode,
                data=birth_wt_df_na_rm,
                ginverse = list(Code=Ainv), # 
                prior = prior,
                pr=TRUE,#saves posterior dist for random effects i.e. what we want 
                nitt=80000,burnin=10000, thin = 50)

plot(birth_wt_model_ROHONLY)


summary(birth_wt_model_ROHONLY)

setwd("/Users/ahewett1/Documents")
save(birth_wt_model_ROHONLY, birth_wt_df_na_rm, file="Rum_deer/ID_Chr_MolEcol_revised/birth_wt_model_output_ROH_ONLY.RData")

