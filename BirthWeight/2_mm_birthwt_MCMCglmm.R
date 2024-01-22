library(tidyverse)
library("MCMCglmm")
library(data.table)
#library(MasterBayes)
library(pedigree)

setwd("H:/")

birth_wt_df=read.table("PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_df.txt", sep=",", header=T)
birth_wt_df_na_rm=birth_wt_df%>%
  na.omit()

ped=read.table("PhD_4th_yr/Inbreeding_depression_models/birth_weight/pedigree.txt", sep=",", header=T)

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

birth_wt_df_na_rm$Code <- as.factor(birth_wt_df_na_rm$Code)
birth_wt_df_na_rm$MumCode <- as.factor(birth_wt_df_na_rm$MumCode)
birth_wt_df_na_rm$MotherStatus <- as.factor(birth_wt_df_na_rm$MotherStatus)
birth_wt_df_na_rm$BirthYear <- as.factor(birth_wt_df_na_rm$BirthYear)
birth_wt_df_na_rm$Sex <- as.factor(birth_wt_df_na_rm$Sex)



k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#### RUNNING MULTI-MEMBERSHIP MODEL ####
birth_wt_model<-MCMCglmm(CaptureWt~1 + Sex + AgeHrs + MotherStatus + mum_age+mum_age_sq+Day_seq+FROHsum, 
                random= ~ Code +
                  idv(FROH_chr1+FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                        FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                        FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                        FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                        FROH_chr33)+BirthYear +MumCode,
                data=birth_wt_df_na_rm,
                ginverse = list(Code=Ainv), # code is recognised as the pedigree 
                prior = prior,
                pr=TRUE,#saves posterior dist for random effects i.e. what we want 
                nitt=100000,burnin=10000, thin = 50)

plot(birth_wt_model)


summary(birth_wt_model)


save(birth_wt_model, file="PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_model_output24072023.RData")

