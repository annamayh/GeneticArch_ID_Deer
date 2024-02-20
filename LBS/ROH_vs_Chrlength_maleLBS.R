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


load("Rum_deer/ID_Chr_MolEcol_revised/male_LBS_model_output_final_ROH_ONLY.RData")

var_LBS=sqrt(sd(male_LBS_df$LBS)) ##variance in LBS

#getting standard devaition in ROH per chr 
sd_mLBS=male_LBS_df%>%
  select(40:72)%>% ## getting only the ROH per chromosome 
  apply(2, sd)%>%
  as_tibble()%>%
  rename(SD_frohchr=value)%>%
  rownames_to_column(var = "CHR")



## solutions for hurdle
FROH_sum_sol_zi=summary(male_LBS_model_ROH_ONLY)$solutions[4,1]
FROH_sum_hu_upr=summary(male_LBS_model_ROH_ONLY)$solutions[4,2]
FROH_sum_hu_lwr=summary(male_LBS_model_ROH_ONLY)$solutions[4,3]

sols_hu<-as.data.frame(male_LBS_model_ROH_ONLY$Sol)%>%dplyr::select(matches("2):ROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_zi)) ## adding FROHsum to chrFROH values

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

size_var_explainedzi=effect_v_sizezi%>%
  inner_join(sd_mLBS)%>%
  mutate(var_exp=(((solution*SD_frohchr)^2)))


effect_chr_maleLBS_zi=ggplot(size_var_explainedzi, aes(x=length_Mb, y=var_exp, label=CHR)) +
  geom_point() +
  geom_smooth(method=lm, color="mediumseagreen") +
  geom_text(nudge_x = 3, size=5)+
  theme_classic()+
  theme(text = element_text(size = 16))+
  labs(x="Chromosome length (Mb)", y="Variance explained by ROHs in zero inflation", tag="D", title = "Male LBS")+
  stat_cor(method = "pearson",size=6)+
  scale_x_continuous(breaks=seq(0, 125, 25))


effect_chr_maleLBS_zi



## same for poisson process

FROH_sum_sol_pois=summary(male_LBS_model_ROH_ONLY)$solutions[3,1]
FROH_sum_pois_upr=summary(male_LBS_model_ROH_ONLY)$solutions[3,2]
FROH_sum_pois_lwr=summary(male_LBS_model_ROH_ONLY)$solutions[3,3]

sols_pois<-as.data.frame(male_LBS_model_ROH_ONLY$Sol)%>%dplyr::select(matches("1):ROH_"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_pois)) ## adding FROHsum to chrFROH values

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

size_var_explained_pois=effect_v_size_pois%>%
  inner_join(sd_mLBS)%>%
  mutate(var_exp=(((solution*SD_frohchr)^2)))


effect_chr_maleLBS_pois=ggplot(size_var_explained_pois, aes(x=length_Mb, y=var_exp, label=CHR)) +
  geom_point() +
  geom_smooth(method=lm, color="chocolate3") +
  geom_text(nudge_x = 3, size=5)+
  theme_classic()+
  theme(text = element_text(size = 16))+
  labs(x="Chromosome length (Mb)", y="Variance explained by ROHs in Poisson")+
  stat_cor(method = "pearson",size=6)+
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
       file = "Rum_deer/ID_Chr_MolEcol_revised/chrlength_Vs_ROH_main_fig.png",
       width = 16,
       height = 13)

