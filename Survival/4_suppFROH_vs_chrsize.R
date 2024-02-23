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
load("PhD/inbreeding_dep_manuscript/model_outputs/A_juvenile_survival_model_output.RData")

sd_js=juvenile_surv_df_na_rm%>%
  select(8:40)%>% ## getting only the ROH per chromosome 79
  apply(2, sd)%>%
  as_tibble()%>%
  rename(SD_frohchr=value)%>%
  rownames_to_column(var = "CHR")

summary(juvenile_surv_model)
ROH_sum_sol=summary(juvenile_surv_model)$solutions[7,1]

rands=as.data.frame(juvenile_surv_model$VCV)

sols_full<-as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("ROH_"))## taking out sols with FROH included

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 

Random_table<-tibble(sols,row.names=names)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "ROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
FROH_sols$CHR<-as.factor(FROH_sols$CHR)


effect_v_size=FROH_sols%>%inner_join(deermap)%>%
  select(CHR, length_Mb, solution)


sd_fixed=juvenile_surv_df_na_rm%>%
  select(c(7, 43:44, 41))%>% 
  apply(2, sd)%>%
  as_tibble()


mum_age=summary(juvenile_surv_model)$solutions[8,1]
mum_age_sq=summary(juvenile_surv_model)$solutions[9,1]
day_seq=summary(juvenile_surv_model)$solutions[10,1]

  
fixed_effects_var=((mum_age*sd_fixed[1,1])^2)+  #for every fixed effect do (beta * sd) which is equivalent to beta^2 * var 
  ((mum_age_sq*sd_fixed[2,1])^2)+ 
  ((day_seq*sd_fixed[3,1])^2)+
  ((summary(juvenile_surv_model)$solutions[1,1])^2)+
  ((ROH_sum_sol*sd_fixed[4,1])^2)

## total variance in survival (including fixed and random effects)
latent_scale_total_var=as.numeric(mean(rowSums(juvenile_surv_model$VCV))+fixed_effects_var)
  

size_var_explained=effect_v_size%>%
  inner_join(sd_js)%>%
  mutate(beta_sd=((solution*SD_frohchr)^2))%>%
  mutate(total_var=(sum(beta_sd)+latent_scale_total_var))%>%
  mutate(var_exp=((solution*SD_frohchr)^2)/total_var)



effect_chr_js_main=ggplot(size_var_explained, aes(x=log(length_Mb), y=log(var_exp), label=CHR)) +
  geom_point(size=13, alpha=0.85) +
  geom_smooth(method=lm, color="red") +
  geom_text(size=7, color = "white")+
  theme_classic()+
  theme(text = element_text(size = 20))+
  labs(x="Log(Chromosome length (Mb))", y="Log(Variance explained by FROH)", title="Juvenile Survival", tag="B")+
  stat_cor(method = "pearson",size=6,label.y.npc = 0.05, label.x.npc = 0.6)#+
  #ylim(0, 0.065)


effect_chr_js_main

# ggsave(effect_chr_js_main,
#        file = "Rum_deer/ID_Chr_MolEcol_revised/juv_surv_var_explained.png",
#        width = 7,
#        height = 5)

sum(size_var_explained$var_exp)
