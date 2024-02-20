

## MALE LBS ###
library(tidyverse)
library(MCMCglmm)


setwd("/Volumes/Seagate Portable Drive")
male_LBS_df=read.table("PhD/Chapter_4_inbr_dep/Male_LBS_df.txt", sep=",", header=T)#692 ids

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

male_LBS_df=male_LBS_df%>%
  left_join(MB_perLG)



male_LBS_df$Code=as.factor(male_LBS_df$Code)
male_LBS_df$BirthYear=as.factor(male_LBS_df$BirthYear)
male_LBS_df$MumCode=as.factor(male_LBS_df$MumCode)

# 
# 


k<-1000
prior.6 = list(R = list(V = diag(2), nu=2, fix = 2),
               G = list(G1 = list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                        G2 = list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                        G3 = list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                        G4 = list(V=1,nu=1,alpha.mu=0,alpha.V=k)))#G3 referring to chrFROH





#################################################################################
# Run model #
#################################################################################
male_LBS_model_ROH_ONLY=MCMCglmm(LBS~ trait+trait:ROHsum-1,
                        
                        random=~idh(trait):MumCode+idh(trait):BirthYear + 
                          idv(at.level(trait,1):(ROH_chr1+ ROH_chr2+ROH_chr3+ROH_chr4+ROH_chr5+ROH_chr6+ROH_chr7+ROH_chr8+
                                ROH_chr9+ROH_chr10+ROH_chr11+ROH_chr12+ROH_chr13+ROH_chr14+ROH_chr15+ROH_chr16+
                                ROH_chr17+ROH_chr18+ROH_chr19+ROH_chr20+ROH_chr21+ROH_chr22+ROH_chr23+ROH_chr24+
                                ROH_chr25+ROH_chr26+ROH_chr27+ROH_chr28+ROH_chr29+ROH_chr30+ROH_chr31+ROH_chr32+
                                ROH_chr33))+
                          idv(at.level(trait,2):(ROH_chr1+ ROH_chr2+ROH_chr3+ROH_chr4+ROH_chr5+ROH_chr6+ROH_chr7+ROH_chr8+
                                                   ROH_chr9+ROH_chr10+ROH_chr11+ROH_chr12+ROH_chr13+ROH_chr14+ROH_chr15+ROH_chr16+
                                                   ROH_chr17+ROH_chr18+ROH_chr19+ROH_chr20+ROH_chr21+ROH_chr22+ROH_chr23+ROH_chr24+
                                                   ROH_chr25+ROH_chr26+ROH_chr27+ROH_chr28+ROH_chr29+ROH_chr30+ROH_chr31+ROH_chr32+
                                                   ROH_chr33)), 
                        rcov = ~ idh(trait):units,
                        data = male_LBS_df, 
                        family = "zipoisson", 
                        prior = prior.6,
                        #verbose = FALSE, 
                        pr = TRUE, 
                        pl = TRUE,
                        nitt=1000000, burnin=250000, thin = 500
)


#first fixed term = intercept of the Poisson process 
#second fixed term = intercept of the zero-inflation  
#remaining terms are trait (Poi & Zi) specific contrasts from the intercept.

## trace plots ##
plot(male_LBS_model_ROH_ONLY)


## model output ##
summary(male_LBS_model_ROH_ONLY)


setwd("/Users/ahewett1/Documents")
save(male_LBS_model_ROH_ONLY,male_LBS_df, file="Rum_deer/ID_Chr_MolEcol_revised/male_LBS_model_output_final_ROH_ONLY.RData")
