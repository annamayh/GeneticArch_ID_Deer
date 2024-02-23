library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)

setwd("/Users/ahewett1/Documents")
deermap <- read.csv("Rum_deer/ID_Chr_MolEcol_revised/mCerEla1.1_assembly.csv", header = T, stringsAsFactors = F)%>%
  select(Chromosome.name, Seq.length)%>%
  rename(CHR=Chromosome.name)%>%
  mutate(length_Mb=Seq.length/1000000)%>%
  filter(!CHR %in% c("Un","X"))



setwd("/Volumes/Seagate Portable Drive")
load("PhD/inbreeding_dep_manuscript/model_outputs/Female_LBS_model_output_final.RData")


#getting standard devaition in ROH per chr 
sd_feLBS=female_LBS_df%>%
  select(4:36)%>% ## getting only the ROH per chromosome 
  apply(2, sd)%>%
  as_tibble()%>%
  rename(SD_frohchr=value)%>%
  rownames_to_column(var = "CHR")

summary(female_LBS_model)


## solutions for hurdle
FROH_sum_sol_hu=summary(female_LBS_model)$solutions[4,1]

sols_hu<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))

names <- sols_hu %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_hu,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_hu,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_hu,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_hu<-Random_table%>%filter(model_variable %like% "ROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
FROH_sols_hu$CHR<-as.factor(FROH_sols_hu$CHR)


effect_v_sizehu=FROH_sols_hu%>%inner_join(deermap)%>%
  select(CHR, length_Mb, solution, CI_upper, CI_lower)



sd_ROH=female_LBS_df%>%
  select(37)%>% 
  apply(2, sd)


fixed_effects_var=
  ((FROH_sum_sol_hu*sd_ROH)^2)

## total variance in survival (including fixed and random effects)
latent_scale_total_var=as.numeric(mean(rowSums(female_LBS_model$VCV))+fixed_effects_var)

size_var_explainedhu=effect_v_sizehu%>%
  inner_join(sd_feLBS)%>%
  mutate(var_exp=(solution*SD_frohchr)^2/latent_scale_total_var)


effect_chr_femaleLBS_hu=ggplot(size_var_explainedhu, aes(x=log(length_Mb), y=log(var_exp), label=CHR)) +
  geom_point(size=10, alpha=0.85) +
  geom_smooth(method=lm, color="mediumseagreen") +
  geom_text(size=6, color = "white")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  labs(x="Log(Chromosome length (Mb))", y="Log(Variance explained by FROH in hurdle)", tag="C", title = "Female LBS")+
  stat_cor(method = "pearson",size=6, label.y.npc = 0.05, label.x.npc = 0.35)+
  scale_x_continuous(breaks=seq(0, 125, 25))


effect_chr_femaleLBS_hu



## same for poisson process

FROH_sum_sol_pois=summary(female_LBS_model)$solutions[3,1]
FROH_sum_pois_upr=summary(female_LBS_model)$solutions[3,2]
FROH_sum_pois_lwr=summary(female_LBS_model)$solutions[3,3]

sols_pois<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))

names <- sols_pois %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_pois,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_pois,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_pois,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_pois<-Random_table%>%filter(model_variable %like% "ROH_c")%>% add_column(CHR = 1:33)#%>% ##filtering all random variables for those including FROH
#add_column(model = "poisson") 
FROH_sols_pois$CHR<-as.factor(FROH_sols_pois$CHR)



effect_v_size_pois=FROH_sols_pois%>%inner_join(deermap)%>%
  select(CHR, length_Mb, solution, CI_upper, CI_lower)

fixed_effects_var_p=
  ((FROH_sum_sol_pois*sd_ROH)^2)

## total variance in survival (including fixed and random effects)
latent_scale_total_var_poi=as.numeric(mean(rowSums(female_LBS_model$VCV))+fixed_effects_var_p)

size_var_explainedhu=effect_v_sizehu%>%
  inner_join(sd_feLBS)%>%
  mutate(var_exp=(solution*SD_frohchr)^2/latent_scale_total_var_poi)

size_var_explained_pois=effect_v_size_pois%>%
  inner_join(sd_feLBS)%>%
  mutate(var_exp=(((solution*SD_frohchr)^2)))


effect_chr_femaleLBS_pois=ggplot(size_var_explained_pois, aes(x=log(length_Mb), y=log(var_exp), label=CHR)) +
  geom_point(size=10, alpha=0.85) +
  geom_smooth(method=lm, color="chocolate3") +
  geom_text(size=6, color = "white")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  labs(x="Log(Chromosome length (Mb))", y="Log(Variance explained by FROH in truncated Poisson)")+
  stat_cor(method = "pearson",size=6, label.y.npc = 0.05, label.x.npc = 0.35)+
  scale_x_continuous(breaks=seq(0, 125, 25))
    


effect_chr_femaleLBS_pois


ROHeff_v_chr_female_LBS=effect_chr_femaleLBS_hu+effect_chr_femaleLBS_pois
ROHeff_v_chr_female_LBS
