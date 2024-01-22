

## MALE LBS ###
library(tidyverse)
library(MCMCglmm)


setwd("H:/")

male_LBS_df=read.table("PhD_4th_yr/Inbreeding_depression_models/LBS/Male_LBS_df.txt", sep=",", header=T)

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
male_LBS_model=MCMCglmm(LBS~ trait+trait:FROHsum-1,
                        
                        random=~idh(trait):MumCode+idh(trait):BirthYear + 
                          idv(at.level(trait,1):(FROH_chr1+ FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                                FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                                FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                                FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                                FROH_chr33))+
                          idv(at.level(trait,2):(FROH_chr1+ FROH_chr2+FROH_chr3+FROH_chr4+FROH_chr5+FROH_chr6+FROH_chr7+FROH_chr8+
                                                   FROH_chr9+FROH_chr10+FROH_chr11+FROH_chr12+FROH_chr13+FROH_chr14+FROH_chr15+FROH_chr16+
                                                   FROH_chr17+FROH_chr18+FROH_chr19+FROH_chr20+FROH_chr21+FROH_chr22+FROH_chr23+FROH_chr24+
                                                   FROH_chr25+FROH_chr26+FROH_chr27+FROH_chr28+FROH_chr29+FROH_chr30+FROH_chr31+FROH_chr32+
                                                   FROH_chr33)), 
                        rcov = ~ idh(trait):units,
                        data = male_LBS_df, 
                        family = "zipoisson", 
                        prior = prior.6,
                        #verbose = FALSE, 
                        pr = TRUE, 
                        pl = TRUE,
                        nitt=1000000, burnin=250000, thin = 500
)


beep(sound = 4)

#first fixed term = intercept of the Poisson process 
#second fixed term = intercept of the zero-inflation  
#remaining terms are trait (Poi & Zi) specific contrasts from the intercept.

## trace plots ##
plot(male_LBS_model)



diag(autocorr(male_LBS_model$VCV)[2,,])


## model output ##
summary(male_LBS_model)



save(male_LBS_model,male_LBS_df, file="PhD_4th_yr/Inbreeding_depression_models/LBS/male_LBS_model_output_final.RData")
