
library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)

setwd("/Users/ahewett1/Documents")
# deer_pc <- read.csv("Rum_deer/ID_Chr_MolEcol_revised/CHR_vs_PC_genes.csv", header = T, stringsAsFactors = F)
# deer_pc$CHR=as.factor(deer_pc$CHR)

deermap <- read.csv("Rum_deer/ID_Chr_MolEcol_revised/mCerEla1.1_assembly.csv", header = T, stringsAsFactors = F)%>%
  select(Chromosome.name, Seq.length)%>%
  rename(CHR=Chromosome.name)%>%
  mutate(length_Mb=Seq.length/1000000)%>%
  filter(!CHR %in% c("Un","X"))

## get sd of FROHchr 
setwd("/Volumes/Seagate Portable Drive")

birth_wt_df_sd=read.table("PhD/Chapter_4_inbr_dep/birth_wt_df.txt", sep=",", header=T)%>%
  na.omit()%>%
  dplyr::select(matches("FROH_chr"))
  
sd_birth_wt_df=apply(birth_wt_df_sd,2, sd)%>%
  as_tibble()%>%
  rename(SD_frohchr=value)%>%
  rownames_to_column(var = "CHR")
  
birth_wt_var=read.table("PhD/Chapter_4_inbr_dep/birth_wt_df.txt", sep=",", header=T)%>%
  na.omit()
  
var_bw=sqrt(sd(birth_wt_var$CaptureWt))

## Read in birth weight output ##
load(file="PhD/inbreeding_dep_manuscript/model_outputs/birth_wt_model_output24072023.RData")


summary_table=summary(birth_wt_model)$solutions
summary_table

FROH_sum_sol=summary(birth_wt_model)$solutions[11,1]
FROHsumest=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(matches("FROHs"))
FROH_sum_hu_upr=quantile(FROHsumest$FROHsum, prob=c(0.975)) 
FROH_sum_hu_lwr=quantile(FROHsumest$FROHsum, prob=c(0.025)) 


rands=as.data.frame(birth_wt_model$VCV)

sols_full<-as.data.frame(birth_wt_model$Sol)%>%dplyr::select(matches("FROH"))%>% ## taking out sols with FROH included
  dplyr::mutate(across(2:34, ~.x + FROH_sum_sol)) ## adding FROHsum to chrFROH values

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
CI_upper<-apply(sols_full,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
CI_lower<-apply(sols_full,2,quantile,probs = c(0.025)) #gets lower CI for all solutions

Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
FROH_sols$CHR<-as.factor(FROH_sols$CHR)


effect_v_size=FROH_sols%>%inner_join(deermap)%>%
  select(CHR, length_Mb, solution, CI_upper, CI_lower)

size_var_explained=effect_v_size%>%
  inner_join(sd_birth_wt_df)%>%
  mutate(var_exp=(((solution*SD_frohchr)^2)/var_bw))
  


effect_chr_bw_main=ggplot(size_var_explained, aes(x=length_Mb, y=var_exp, label=CHR)) +
  geom_point() +
  geom_smooth(method=lm, color="red") +
  geom_text(nudge_x = 3, size=5)+
  theme_classic()+
  theme(text = element_text(size = 20))+
  labs(x="Chromosome length (Mb)", y="ID variance explained", title="Birth weight", tag="A")+
  stat_cor(method = "pearson",size=6)
  
  effect_chr_bw_main
  
  # effect_chr_bw_errors=ggplot(effect_v_size, aes(x=length_Mb, y=solution, label=CHR, ymin=CI_lower, ymax=CI_upper)) +
  #   geom_pointrange() +
  #   geom_smooth(method=lm, color="red") +
  #   geom_text(nudge_x = 3, size=5)+
  #   theme_classic()+
  #   theme(text = element_text(size = 22))+
  #   labs(x="Chromosome length (Mb)", y="Slope estimate", title="Birth weight", tag="A")+
  #   stat_cor(method = "pearson",label.x = 110, size=6)
  # 
  # effect_chr_bw_errors

######################################################################################################

## with survival ##
  
  
surv_df_sd=read.table("PhD/Chapter_4_inbr_dep/AA_juvenile_survival_df.txt", sep=",", header=T)%>%
    na.omit()%>%
    dplyr::select(matches("FROH_chr"))
  
sd_surv_df=apply(surv_df_sd,2, sd)%>%
    as_tibble()%>%
    rename(SD_frohchr=value)%>%
    rownames_to_column(var = "CHR")
  

load(file="PhD/inbreeding_dep_manuscript/model_outputs/A_juvenile_survival_model_output.RData")
  
  summary(juvenile_surv_model)
  
  summary_table=summary(juvenile_surv_model)$solutions
  summary_table
  
  FROH_sum_sol=summary(juvenile_surv_model)$solutions[7,1]
  FROH_sum_upr=summary(juvenile_surv_model)$solutions[7,2]
  FROH_sum_lwr=summary(juvenile_surv_model)$solutions[7,3]
  
  
  mean(rowSums(juvenile_surv_model$VCV))
  
  
  sols_full<-as.data.frame(juvenile_surv_model$Sol)%>%dplyr::select(matches("FROH_"))%>% ## taking out sols with FROH included
    dplyr::mutate(across(1:33, ~.x + FROH_sum_sol)) ## adding FROHsum to chrFROH values
  
  names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
  sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(sols_full,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(sols_full,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)
  
  names(Random_table)[1]<-"solution"
  names(Random_table)[2]<-"model_variable"
  
  FROH_sols<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
  
  FROH_sols$CHR<-as.factor(FROH_sols$CHR)


  
  effect_v_size=FROH_sols%>%inner_join(deermap)%>%
    select(CHR, length_Mb, solution, CI_upper, CI_lower)
  
  
  size_var_explained=effect_v_size%>%
    inner_join(sd_surv_df)%>%
    mutate(var_exp=((solution*SD_frohchr)^2))
  
  
  effect_chr_surv_main=ggplot(size_var_explained, aes(x=length_Mb, y=var_exp, label=CHR)) +
    geom_point() +
    geom_smooth(method=lm, color="red") +
    geom_text(nudge_x = 3, size=5)+
    theme_classic()+
    theme(text = element_text(size = 22))+
    labs(x="Chromosome length (Mb)", y="ID variance explained", title="Juvenile survival", tag="B")+
    stat_cor(method = "pearson",size=6)
  
  effect_chr_surv_main
  
  
  
  # effect_chr_surv_errors=ggplot(effect_v_size, aes(x=length_Mb, y=solution, label=CHR, ymin=CI_lower, ymax=CI_upper)) +
  #   geom_pointrange() +
  #   geom_smooth(method=lm, color="red") +
  #   geom_text(nudge_x = 3, size=5)+
  #   theme_classic()+
  #   theme(text = element_text(size = 22))+
  #   labs(x="Chromosome length (Mb)", y="Slope estimate", title="Juvenile survival", tag="B")+
  #   stat_cor(method = "pearson",label.x = 110, size=6)
  # 
  # effect_chr_surv_errors




  #############################################################################################################

  ## female LBS ##
  
  
  femaleLBS_df_sd=read.table("PhD/Chapter_4_inbr_dep/female_LBS_df.txt", sep=",", header=T)%>%
    na.omit()%>%
    dplyr::select(matches("FROH_chr"))
  
  sd_femaleLBS=apply(femaleLBS_df_sd,2, sd)%>%
    as_tibble()%>%
    rename(SD_frohchr=value)%>%
    rownames_to_column(var = "CHR")
  
  
  load(file="PhD/inbreeding_dep_manuscript/model_outputs/Female_LBS_model_output_final.RData")
  
  
  summary(female_LBS_model)
  #plot(male_LBS_model)
  
  FROH_sum_sol_hu=summary(female_LBS_model)$solutions[4,1]
  FROH_sum_hu_upr=summary(female_LBS_model)$solutions[4,2]
  FROH_sum_hu_lwr=summary(female_LBS_model)$solutions[4,3]
  
  
  sols_hu<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))%>% ## taking out sols with FROH included
    dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_hu)) ## adding FROHsum to chrFROH values
  
  names <- sols_hu %>% names ## gets names of all random variables, 2 = all down row
  sols<-apply(sols_hu,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(sols_hu,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(sols_hu,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)
  
  names(Random_table)[1]<-"solution"
  names(Random_table)[2]<-"model_variable"
  
  FROH_sols_hu<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH

  
  
  ## same plot for poisson process
  
  FROH_sum_sol_pois=summary(female_LBS_model)$solutions[3,1]
  FROH_sum_pois_upr=summary(female_LBS_model)$solutions[3,2]
  FROH_sum_pois_lwr=summary(female_LBS_model)$solutions[3,3]
  
  sols_pois<-as.data.frame(female_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))%>% ## taking out sols with FROH included
    dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_pois)) ## adding FROHsum to chrFROH values
  
  names <- sols_pois %>% names ## gets names of all random variables, 2 = all down row
  sols<-apply(sols_pois,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(sols_pois,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(sols_pois,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)
  
  names(Random_table)[1]<-"solution"
  names(Random_table)[2]<-"model_variable"
  
  FROH_sols_pois<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33)#%>% ##filtering all random variables for those including FROH
  #add_column(model = "poisson") 
  
  FROH_sols_hu$CHR<-as.factor(FROH_sols_hu$CHR)
  

  effect_v_size_hu=FROH_sols_hu%>%inner_join(deermap)%>%
    select(CHR, length_Mb,solution,CI_upper, CI_lower)
  
  
  size_var_explained_hu=effect_v_size_hu%>%
    inner_join(sd_femaleLBS)%>%
    mutate(var_exp=((solution*SD_frohchr)^2))
  
# zi_error=ggplot(effect_v_size_hu, aes(x=length_Mb, y=solution, label=CHR,ymin=CI_lower, ymax=CI_upper)) +
#     geom_pointrange() +
#     geom_smooth(method=lm, color="mediumseagreen") +
#     geom_text(nudge_x = 5,size=5)+
#     theme_classic()+
#     theme(text = element_text(size = 22), axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
#     labs(x="Chromosome length (Mb)", y="Slope estimate for hurdle", tag="C", title = "Female LBS")+
#     stat_cor(method = "pearson",size=5)+
#   scale_x_continuous(breaks=seq(0, 125, 25))
# 
#   
  FROH_sols_pois$CHR<-as.factor(FROH_sols_pois$CHR)
  
  effect_v_sizepi=FROH_sols_pois%>%inner_join(deermap)%>%
    select(CHR, length_Mb,solution,CI_upper, CI_lower)

  
  size_var_explained_pi=effect_v_sizepi%>%
    inner_join(sd_femaleLBS)%>%
    mutate(var_exp=((solution*SD_frohchr)^2))  
  
    # 
  # pi_error=ggplot(effect_v_sizepi, aes(x=length_Mb, y=solution, label=CHR,ymin=CI_lower, ymax=CI_upper)) +
  #   geom_pointrange() +
  #   geom_smooth(method=lm, color="chocolate3") +
  #   geom_text(nudge_x = 5,size=5)+
  #   theme_classic()+
  #   theme(text = element_text(size = 22),  axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
  #   labs(x="Chromosome length (Mb)", y="Slope estimate for Poisson process")+
  #   stat_cor(method = "pearson",size=5)+
  #   scale_x_continuous(breaks=seq(0, 125, 25))
  # 
  # 
  # 
  # eff_v_chr_female_LBS_error=zi_error+pi_error
  # eff_v_chr_female_LBS_error
  
  
  zi=ggplot(size_var_explained_hu, aes(x=length_Mb, y=var_exp, label=CHR)) +
    geom_point() +
    geom_smooth(method=lm, color="mediumseagreen") +
    geom_text(nudge_x = 5,size=5)+
    theme_classic()+
    theme(text = element_text(size = 22), axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
    labs(x="Chromosome length (Mb)", y="ID variance explained hurdle", tag="C", title = "Female LBS")+
    stat_cor(method = "pearson",size=5, label.y = 165)+
    scale_x_continuous(breaks=seq(0, 125, 25))
  
  
  
  pi=ggplot(size_var_explained_pi, aes(x=length_Mb, y=var_exp, label=CHR)) +
    geom_point() +
    geom_smooth(method=lm, color="chocolate3") +
    geom_text(nudge_x = 5,size=5)+
    theme_classic()+
    theme(text = element_text(size = 22),  axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
    labs(x="Chromosome length (Mb)", y="ID variance explained Poisson")+
    stat_cor(method = "pearson",size=5, label.y = 0.5) +
    scale_x_continuous(breaks=seq(0, 125, 25))
  

  
  
  eff_v_chr_female_LBS=zi+pi
  eff_v_chr_female_LBS

#############################################################################################################
  
  
  
  maleLBS_df_sd=read.table("PhD/Chapter_4_inbr_dep/Male_LBS_df.txt", sep=",", header=T)%>%
    na.omit()%>%
    dplyr::select(matches("FROH_chr"))
  
  sd_maleLBS=apply(maleLBS_df_sd,2, sd)%>%
    as_tibble()%>%
    rename(SD_frohchr=value)%>%
    rownames_to_column(var = "CHR")
  
  setwd("/Volumes/Seagate Portable Drive")
  load(file="PhD/Chapter_4_inbr_dep/male_LBS_model_output_final.RData")
  
  
  summary(male_LBS_model)
  #plot(male_LBS_model)
  
  FROH_sum_sol_zi=summary(male_LBS_model)$solutions[4,1]
  FROH_sum_zi_upr=summary(male_LBS_model)$solutions[4,2]
  FROH_sum_zi_lwr=summary(male_LBS_model)$solutions[4,3]
  
  
  sols_zi<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))%>% ## taking out sols with FROH included
    dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_zi)) ## adding FROHsum to chrFROH values
  
  names <- sols_zi %>% names ## gets names of all random variables, 2 = all down row
  sols<-apply(sols_zi,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(sols_zi,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(sols_zi,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)
  
  names(Random_table)[1]<-"solution"
  names(Random_table)[2]<-"model_variable"
  
  FROH_sols_zi<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
  
  FROH_sum_sol_pois=summary(male_LBS_model)$solutions[3,1]
  FROH_sum_pois_upr=summary(male_LBS_model)$solutions[3,2]
  FROH_sum_pois_lwr=summary(male_LBS_model)$solutions[3,3]
  
  sols_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))%>% ## taking out sols with FROH included
    dplyr::mutate(across(1:33, ~.x + FROH_sum_sol_pois)) ## adding FROHsum to chrFROH values
  
  names <- sols_pois %>% names ## gets names of all random variables, 2 = all down row
  sols<-apply(sols_pois,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(sols_pois,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(sols_pois,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  Random_table<-tibble(sols,row.names=names)%>%add_column(CI_upper)%>%add_column(CI_lower)
  
  names(Random_table)[1]<-"solution"
  names(Random_table)[2]<-"model_variable"
  
  FROH_sols_pois<-Random_table%>%filter(model_variable %like% "FROH_c")%>% add_column(CHR = 1:33)#%>% ##filtering all random variables for those including FROH
  #add_column(model = "poisson") 

  
  FROH_sols_zi$CHR<-as.factor(FROH_sols_zi$CHR)
  
  
  effect_v_size_zi=FROH_sols_zi%>%inner_join(deermap)%>%
    select(CHR, length_Mb,solution,CI_upper, CI_lower)
  
  # zi_error_m=ggplot(effect_v_size_zi, aes(x=length_Mb, y=solution, label=CHR,ymin=CI_lower, ymax=CI_upper)) +
  #   geom_pointrange() +
  #   geom_smooth(method=lm, color="mediumseagreen") +
  #   geom_text(nudge_x = 5,size=5)+
  #   theme_classic()+
  #   theme(text = element_text(size = 22), axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
  #   labs(x="Chromosome length (Mb)", y="Slope estimate for zero inflation", tag="D", title = "Male LBS")+
  #   stat_cor(method = "pearson",size=5)+
  #   scale_x_continuous(breaks=seq(0, 125, 25))
  # 
  # 
  FROH_sols_pois$CHR<-as.factor(FROH_sols_pois$CHR)
  # 
  effect_v_sizepi_m=FROH_sols_pois%>%inner_join(deermap)%>%
    select(CHR, length_Mb,solution,CI_upper, CI_lower)
  # 
  # pi_error_m=ggplot(effect_v_sizepi_m, aes(x=length_Mb, y=solution, label=CHR,ymin=CI_lower, ymax=CI_upper)) +
  #   geom_pointrange() +
  #   geom_smooth(method=lm, color="chocolate3") +
  #   geom_text(nudge_x = 5,size=5)+
  #   theme_classic()+
  #   theme(text = element_text(size = 22),  axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
  #   labs(x="Chromosome length (Mb)", y="Slope estimate for Poisson process")+
  #   stat_cor(method = "pearson",size=5)+
  #   scale_x_continuous(breaks=seq(0, 125, 25))
  # 
  # 
  
  # eff_v_chr_male_LBS_error=zi_error_m+pi_error_m
  # eff_v_chr_male_LBS_error
 
  
  
  size_var_explained_zi=effect_v_size_zi%>%
    inner_join(sd_maleLBS)%>%
    mutate(var_exp=((solution*SD_frohchr)^2))
  
  
  zi_m=ggplot(size_var_explained_zi, aes(x=length_Mb, y=var_exp, label=CHR)) +
    geom_point() +
    geom_smooth(method=lm, color="mediumseagreen") +
    geom_text(nudge_x = 5,size=5)+
    theme_classic()+
    theme(text = element_text(size = 22), axis.title.x = element_text(size=12),axis.text.x = element_text(size=12))+
    labs(x="Chromosome length (Mb)", y="Slope estimate for zero inflation", tag="D", title = "Male LBS")+
    stat_cor(method = "pearson",size=5, label.y = 150)+
    scale_x_continuous(breaks=seq(0, 125, 25))
  
  
  
  size_var_explained_pi=effect_v_sizepi%>%
    inner_join(sd_maleLBS)%>%
    mutate(var_exp=((solution*SD_frohchr)^2))
  
  
  pi_m=ggplot(size_var_explained_pi, aes(x=length_Mb, y=var_exp, label=CHR)) +
    geom_point() +
    geom_smooth(method=lm, color="chocolate3") +
    geom_text(nudge_x = 5,size=5)+
    theme_classic()+
    theme(text = element_text(size = 22),  axis.title.x = element_text(size=12), axis.text.x = element_text(size=12))+
    labs(x="Chromosome length (Mb)", y="Slope estimate for Poisson process")+
    stat_cor(method = "pearson",size=5, label.y = 0.6) +
    scale_x_continuous(breaks=seq(0, 125, 25))
  
  
  
  eff_v_chr_male_LBS=zi_m+pi_m
  eff_v_chr_male_LBS
  

  
bottom_main=zi+pi+zi_m+pi_m+plot_layout(nrow=1)
bottom_main

main_fig=(effect_chr_bw_main+effect_chr_surv_main)/bottom_main+plot_layout(nrow=2)
main_fig



# bottom_supp=zi_error+pi_error+zi_error_m+pi_error_m+plot_layout(nrow=1)
# bottom_supp
# 
# supp_fig=(effect_chr_bw_errors+effect_chr_surv_errors)/bottom_supp+plot_layout(nrow=2)
# supp_fig


setwd("/Users/ahewett1/Documents")

ggsave(main_fig,
       file = "Rum_deer/ID_Chr_MolEcol_revised/chr_size_slope_main_fig.png",
       width = 16,
       height = 13)

# ggsave(supp_fig,
#        file = "Rum_deer/ID_Chr_MolEcol_revised/SUPP_chr_size_slop_fig_errors.png",
#        width = 16,
#        height = 13)


#############################################################################################################
# pc_vs_length=deer_pc%>%
#   inner_join(deermap)
# 
# ggplot(pc_vs_length, aes(x=Num_PC_genes, y=length_Mb, label=CHR)) +
#   geom_point() +
#   geom_smooth(method=lm, color="red") +
#   geom_text(nudge_x = 30)+
#   theme_classic()+
#   theme(text = element_text(size = 18))+
#   labs(x="Number of protein coding genes", y="Chromosome size")+
#   stat_cor(method = "pearson",label.x = 110)
