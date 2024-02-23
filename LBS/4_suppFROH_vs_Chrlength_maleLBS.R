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
load("PhD/inbreeding_dep_manuscript/model_outputs/male_LBS_model_output_final.RData")

#getting standard devaition in ROH per chr 
sd_mLBS=male_LBS_df%>%
  select(4:36)%>% ## getting only the ROH per chromosome 
  apply(2, sd)%>%
  as_tibble()%>%
  rename(SD_frohchr=value)%>%
  rownames_to_column(var = "CHR")

summary(male_LBS_model)

## solutions for hurdle
FROH_sum_sol_zi=summary(male_LBS_model)$solutions[4,1]


sols_hu<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))
names <- sols_hu %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_hu,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_hu,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_hu,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols_zi<-Random_table%>%filter(model_variable %like% "ROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
FROH_sols_zi$CHR<-as.factor(FROH_sols_zi$CHR)


effect_v_sizezi=FROH_sols_zi%>%inner_join(deermap)%>%
  select(CHR, length_Mb, solution, CI_upper, CI_lower)


sd_ROH=male_LBS_df%>%
  select(37)%>% 
  apply(2, sd)


fixed_effects_var=
  ((FROH_sum_sol_zi*sd_ROH)^2)

## total variance in survival (including fixed and random effects)
latent_scale_total_var=as.numeric(mean(rowSums(male_LBS_model$VCV))+fixed_effects_var)


size_var_explainedzi=effect_v_sizezi%>%
  inner_join(sd_mLBS)%>%
  mutate(var_exp=(((solution*SD_frohchr)^2)/latent_scale_total_var))


effect_chr_maleLBS_zi=ggplot(size_var_explainedzi, aes(x=log(length_Mb), y=log(var_exp), label=CHR)) +
  geom_point(size=10, alpha=0.85) +
  geom_smooth(method=lm, color="mediumseagreen") +
  geom_text(size=6, color = "white")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  labs(x="Log(Chromosome length (Mb))", y="Log(Variance explained by FROH in zero inflation)", tag="D", title = "Male LBS")+
  stat_cor(method = "pearson",size=6,label.y.npc = 0.05, label.x.npc = 0.35)+
  scale_x_continuous(breaks=seq(0, 125, 25))


effect_chr_maleLBS_zi



## same for poisson process

FROH_sum_sol_pois=summary(male_LBS_model)$solutions[3,1]

sols_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))

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



fixed_effects_var_poi=
  ((FROH_sum_sol_pois*sd_ROH)^2)

## total variance in survival (including fixed and random effects)
latent_scale_total_var=as.numeric(mean(rowSums(male_LBS_model$VCV))+fixed_effects_var_poi)


size_var_explained_pois=effect_v_size_pois%>%
  inner_join(sd_mLBS)%>%
  mutate(var_exp=(((solution*SD_frohchr)^2)/latent_scale_total_var))


effect_chr_maleLBS_pois=ggplot(size_var_explained_pois, aes(x=log(length_Mb), y=log(var_exp), label=CHR)) +
  geom_point(size=10, alpha=0.85) +
  geom_smooth(method=lm, color="chocolate3") +
  geom_text(size=6, color = "white")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  labs(x="Log(Chromosome length (Mb))", y="Log(Variance explained by FROH in Poisson)")+
  stat_cor(method = "pearson",size=6, label.y.npc = 0.05, label.x.npc = 0.35)+
  scale_x_continuous(breaks=seq(0, 125, 25))


effect_chr_maleLBS_pois




ROHeff_v_chr_male_LBS=effect_chr_maleLBS_zi+effect_chr_maleLBS_pois
ROHeff_v_chr_male_LBS






###### all others need to br run through first

bottom_main=effect_chr_femaleLBS_hu+effect_chr_femaleLBS_pois+effect_chr_maleLBS_zi+effect_chr_maleLBS_pois+plot_layout(nrow=1)
bottom_main

main_fig=(effect_chr_bw_main+effect_chr_js_main)/bottom_main+plot_layout(nrow=2)
main_fig



setwd("/Users/ahewett1/Documents")

ggsave(main_fig,
       file = "Rum_deer/ID_Chr_MolEcol_revised/supp_chrlength_Vs_FROH.png",
       width = 20,
       height = 13)

