library(tidyverse)
library("MCMCglmm")
library(data.table)

setwd("/Volumes/Seagate Portable Drive")
juvenile_surv_df=read.table("PhD/Chapter_4_inbr_dep/AA_juvenile_survival_df.txt", sep=",", header=T)

juvenile_surv_df_na_rm=juvenile_surv_df%>%
  select(-BirthWt)%>%
  na.omit() #remove all NAs quite a few for BW and mother stat


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

juvenile_surv_df_na_rm=juvenile_surv_df_na_rm%>%
  left_join(MB_perLG)


## change some variables to factors 

juvenile_surv_df_na_rm$MotherStatus[juvenile_surv_df_na_rm$MotherStatus=='Na\xefve']<-'Naive'


juvenile_surv_df_na_rm$BirthYear=as.factor(juvenile_surv_df_na_rm$BirthYear)
juvenile_surv_df_na_rm$Code=as.factor(juvenile_surv_df_na_rm$Code)
juvenile_surv_df_na_rm$MumCode=as.factor(juvenile_surv_df_na_rm$MumCode)
juvenile_surv_df_na_rm$MotherStatus=as.factor(juvenile_surv_df_na_rm$MotherStatus)
juvenile_surv_df_na_rm$Sex=as.factor(juvenile_surv_df_na_rm$Sex)

k<-100
prior<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k), ## multimemberhsip part
                   G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k)))



juvenile_surv_model_ROH_ONLY<-MCMCglmm(juvenile_survival~1 + Sex + MotherStatus + ROH_sum_div + mum_age+mum_age_sq+Day_seq, 
                 random= ~ 
                   idv(ROH_chr1+ROH_chr2+ROH_chr3+ROH_chr4+ROH_chr5+ROH_chr6+ROH_chr7+ROH_chr8+
                         ROH_chr9+ROH_chr10+ROH_chr11+ROH_chr12+ROH_chr13+ROH_chr14+ROH_chr15+ROH_chr16+
                         ROH_chr17+ROH_chr18+ROH_chr19+ROH_chr20+ROH_chr21+ROH_chr22+ROH_chr23+ROH_chr24+
                         ROH_chr25+ROH_chr26+ROH_chr27+ROH_chr28+ROH_chr29+ROH_chr30+ROH_chr31+ROH_chr32+
                         ROH_chr33)+BirthYear +MumCode,
                 family="threshold",
                 data=juvenile_surv_df_na_rm,
                 prior = prior,
                 pr=TRUE,
                 nitt=300000,burnin=50000, thin = 100
                 )##

# 200k finished in <20 mins 

plot(juvenile_surv_model_ROH_ONLY)

summary(juvenile_surv_model_ROH_ONLY)

setwd("/Users/ahewett1/Documents")
save(juvenile_surv_model_ROH_ONLY, juvenile_surv_df_na_rm, file="Rum_deer/ID_Chr_MolEcol_revised/Juvenile_survival_model_output_ROH_ONLY.RData")


