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


load("Rum_deer/ID_Chr_MolEcol_revised/ROH_ONLY_model_outputs/birth_wt_model_output_ROH_ONLY.RData")

var_bw=(var(birth_wt_df_na_rm$CaptureWt))


sd_bw=birth_wt_df_na_rm%>%
  select(c(47:79))%>% ## getting only the ROH per chromosome 
  apply(2, sd)%>%
  as_tibble()%>%
  rename(SD_frohchr=value)%>%
  rownames_to_column(var = "CHR")

cov=birth_wt_df_na_rm%>%
  select(c(47:79), 81)%>%
  cov()

ROH_sum_sol=summary(birth_wt_model_ROHONLY)$solutions[11,1]

rands=as.data.frame(birth_wt_model_ROHONLY$VCV)

sols_full<-as.data.frame(birth_wt_model_ROHONLY$Sol)%>%select(matches("ROH_"))

names <- sols_full %>% names ## gets names of all random variables, 2 = all down row
sols<-apply(sols_full,2,mean)#gets mean of all solutions i.e. the effect size of random effects 

Random_table<-tibble(sols,row.names=names)

names(Random_table)[1]<-"solution"
names(Random_table)[2]<-"model_variable"

FROH_sols<-Random_table%>%filter(model_variable %like% "ROH_c")%>% add_column(CHR = 1:33) ##filtering all random variables for those including FROH
FROH_sols$CHR<-as.factor(FROH_sols$CHR)


effect_v_size=FROH_sols%>%inner_join(deermap)%>%
  select(CHR, length_Mb, solution)

size_var_explained=effect_v_size%>%
  inner_join(sd_bw)%>%
  mutate(var_exp=((solution*SD_frohchr)^2)/var_bw)


effect_chr_bw_main=ggplot(size_var_explained, aes(x=log(length_Mb), y=log(var_exp), label=CHR)) +
  geom_point(size=13, alpha=0.85) +
  geom_smooth(method=lm, color="red") +
  geom_text(size=7, color = "white")+
  theme_classic()+
  theme(text = element_text(size = 20))+
  labs(x="Log(Chromosome length (Mb))", y="Log(Variance explained by ROHs)", title="Birth weight", tag="A")+
  stat_cor(method = "pearson",size=6,label.y.npc = 0.05, label.x.npc = 0.6)

effect_chr_bw_main


sum(size_var_explained$var_exp)

