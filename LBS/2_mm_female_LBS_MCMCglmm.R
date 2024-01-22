
library(MCMCglmm)
library(data.table)

setwd("H:/")

female_LBS_df=read.table("PhD_4th_yr/Inbreeding_depression_models/LBS/female_LBS_df.txt", sep=",", header=T)#692 ids
female_LBS_df$Code=as.factor(female_LBS_df$Code)
female_LBS_df$BirthYear=as.factor(female_LBS_df$BirthYear)
female_LBS_df$MumCode=as.factor(female_LBS_df$MumCode)


k<-1000
prior.6 = list(R = list(V = diag(2), nu=2, fix = 2),
               G = list(G1 = list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                        G2 = list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                        G3 = list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                        G4 = list(V=1,nu=1,alpha.mu=0,alpha.V=k)))#G3 referring to chrFROH





female_LBS_model=MCMCglmm(LBS~ trait+trait:FROHsum-1,
                          
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
                          data = female_LBS_df, 
                          family = "hupoisson", 
                          prior = prior.6,
                          pr = TRUE, 
                          pl = TRUE,
                          nitt=650000, burnin=150000, thin = 200
)


beep(sound=8)


plot(female_LBS_model)


summary(female_LBS_model)

diag(autocorr(female_LBS_model$VCV)[2,,])
##still a bit of suto in mmcode when thin=200


save(female_LBS_df,female_LBS_model, file="PhD_4th_yr/Inbreeding_depression_models/LBS/Female_LBS_model_output_final.RData")
