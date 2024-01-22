

## jarrods function for converting male LBS to total LBS 

library(tidyverse)
library(data.table)
library(MCMCglmm)
library(patchwork)
library(grid)
library(gridExtra)


setwd("H:/")

load(file="PhD_4th_yr/Inbreeding_depression_models/LBS/male_LBS_model_output_final.RData")
load("PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/plot_dfs_male_LBS.RData") ##load in df made using loops in this script as they take ages!


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


male_LBS_zi=
  ggplot(data=FROH_sols_zi, aes(x=as.factor(CHR), y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate + CI ",title = "Zero inflation", tag = "B") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none",
        text = element_text(size = 18), plot.title = element_text(size=14), 
        axis.text.y = element_text(size=12),axis.title.x = element_text(size=10))+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_zi, yend = FROH_sum_sol_zi, colour = "mediumseagreen", alpha=0.6, linewidth=1)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol_zi, ymin = FROH_sum_zi_upr, ymax = FROH_sum_zi_lwr,
           colour = "mediumseagreen", linewidth = 1, alpha=0.5, size=0.2)+
  expand_limits(x = 35)+
  geom_text(aes(x=34.5, y=FROH_sum_sol_zi, label="*"), colour="mediumseagreen", alpha=0.5)
  

male_LBS_zi


## same plot for poisson process

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

male_LBS_pois=
  ggplot(data=FROH_sols_pois, aes(x=as.factor(CHR), y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange() + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior mean estimate + CI ", title = "Poisson process") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none",
        text = element_text(size = 18), plot.title = element_text(size=14), axis.text.y = element_text(size=12),
        axis.title.y = element_blank(), axis.title.x = element_text(size=10))+
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol_pois, yend = FROH_sum_sol_pois, colour = "chocolate3", alpha=0.6, linewidth=1)+
  ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol_pois, ymin = FROH_sum_pois_upr, ymax = FROH_sum_pois_lwr,
           colour = "chocolate3", linewidth = 1, alpha=0.5, size=0.2)+

  geom_text(aes(x=34.5, y=FROH_sum_sol_pois, label="***"), colour="chocolate3", alpha=0.5)+
  expand_limits(x = 35)

male_LBS_pois



forest=(male_LBS_zi+male_LBS_pois)


# p4 <- ggplot(data.frame(l = male_LBS_zi$labels$y, x = 1, y = 1)) +
#   geom_text(aes(x, y, label = l), angle = 0, size=6) + 
#   theme_void() +
#   coord_cartesian(clip = "off")
# 
# all_forest=(forest/p4) + plot_layout(heights = c(25,1))
# all_forest
# 


###################################################################################
### transforming the data ########################################################
##################################################################################


#### function written by jarrod to tranform the data #### 
# Mu takes two values: the untrsanformed equation output for the poisson part followed by the untrsanformed 
# equation output for the Bernoulli part. V is a 2x2 covariance matrix describing the variance of the 
# random effects you want to average over. Because you fitted idh structures V is a diagonal matrix. 
# You can sum the variances over random terms (mum, year AND residuals) to get V.

normal.zipooisson<-function(mu, V){
  
  int.foo<-function(x, mu, V){(1-plogis(x[2]))*exp(x[1])*mvtnorm::dmvnorm(x, mu, V)}
  
  
  cubature::adaptIntegrate(int.foo, qnorm(0.0001, mu,sqrt(diag(V))), qnorm(0.9999, mu,sqrt(diag(V))), mu=mu, V=V)[[1]]
  
}

##V<-rIW(diag(2), 10) example


summary_table=summary(male_LBS_model)$solutions
summary_table

summary(male_LBS_model)

# Iterations = 250001:999501
# Thinning interval  = 500
# Sample size  = 1500 
# 
# DIC: 1260.806 
# 
# G-structure:  ~idh(trait):MumCode
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.MumCode       0.2632 5.195e-08   0.6113   1500.0
# traitzi_LBS.MumCode    0.4750 9.773e-08   1.7112    688.9
# 
# ~idh(trait):BirthYear
# 
# post.mean  l-95% CI u-95% CI eff.samp
# traitLBS.BirthYear      0.05643 1.469e-11   0.2098     1500
# traitzi_LBS.BirthYear   3.53473 1.289e+00   6.1533     1079

# R-structure:  ~idh(trait):units
# 
# post.mean l-95% CI u-95% CI eff.samp
# traitLBS.units       0.6036   0.2371   0.9747     1500
# traitzi_LBS.units    1.0000   1.0000   1.0000        0



rand_var_zi=0.4750+3.53473+1.0000 ##variance estimatd for radoms effects in zi: mum code, birth year and residual

rand_var_pois=0.2632+0.05643+0.6036 ##variance estimatd for all radoms effects in pois

v=diag(c(rand_var_pois,rand_var_zi),2) ##assuming no covariance 


zi_mean_chr<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("2):FROH_"))%>%## taking out sols with FROH included
  apply(2, mean)

pois_mean_chr<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(matches("1):FROH_"))%>%## taking out sols with FROH included
  apply(2, mean)

#as.numeric(zi_mean_chr[1])

# ibc=0.1
# 
# chr=1

zi_inter_mean=summary_table[2,1]
zi_FROHsum_mean=summary_table[4,1]


pois_inter_mean=summary_table[1,1]
pois_FROHsum_mean=summary_table[3,1]

ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)

LBS_male_full=list()

for (c in 1:length(ibc_qua)){
  
  ibc=ibc_qua[c]
  
  LBS_male=list()
  
  for (chr in 1:33){

    untrans_LBS_zi=(zi_inter_mean+(ibc*(as.numeric(FROH_sols_zi[chr,1])))) 
    untrans_LBS_pois=(pois_inter_mean+(ibc*(as.numeric(FROH_sols_pois[chr,1]))))
    mu=(c(untrans_LBS_pois,untrans_LBS_zi))#  zi and pois effects

    LBS_trans=normal.zipooisson(mu,v) # function written by jarrod taking in to account variance 
    LBS_trans_mat=matrix(c(LBS_trans, paste0(chr), paste0(ibc)), nrow = 1)
    
    LBS_male[[chr]]=LBS_trans_mat
      }

  LBS_male_all_list=do.call(rbind.data.frame,LBS_male)
  LBS_male_full[[c]]=LBS_male_all_list
  
}

beep(sound=1)

LBS_male_all=do.call(rbind.data.frame,LBS_male_full)

names(LBS_male_all)=(c("LBS", "CHR", "ibc"))

LBS_male_all <- sapply(LBS_male_all, as.numeric)%>%as.data.frame()


LBS_pred_independent=ggplot(data=LBS_male_all,aes(x=ibc,y=LBS, group=as.factor(CHR), color=as.factor(CHR)))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding Coefficient", y="LBS", colour="Chromosome", title="Assuming chromosome independence",tag = "C")+
  theme(text = element_text(size = 18), plot.title = element_text(size=13), 
        legend.title = element_text(size=11), 
        legend.text = element_text(size=9))


LBS_pred_independent

##add in real data too?


###########################################################################################################
# for all iterations get mean LBS and 95% CIs 
######################################################################################

sol_inter_pois<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(("traitLBS"))## taking out sols of intercept
sol_inter_zi<-as.data.frame(male_LBS_model$Sol)%>%dplyr::select(("traitzi_LBS"))## taking out sols of intercept

rands_pois=as.data.frame(male_LBS_model$VCV)%>%dplyr::select(matches("traitLBS"))
rands_zi=as.data.frame(male_LBS_model$VCV)%>%dplyr::select(matches("traitzi"))



#rand_var_zi=0.4750+3.53473+1.0000  ##variance estimatd for radoms effects in zi: mum code, birth year and residual
#rand_var_pois=0.2632+0.05643+0.6036 ##variance estimatd for all radoms effects in pois
#v=diag(c(rand_var_pois,rand_var_zi),2) ##assuming no covariance 

ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)
LBS_transformed=list()


for (c in 1:length(ibc_qua)){
  
  ibc=ibc_qua[c]
  LBS_sum=list()
  
  for (row in 1:nrow(sols_pois)){##for every iteration 
    
    rand_var_zi=rands_zi[row,1]+rands_zi[row,2]+rands_zi[row,3] ##variance estimatd for radoms effects in zi: mum code, birth year and residual
    rand_var_pois=rands_pois[row,1]+rands_pois[row,2]+rands_pois[row,3]  ##variance estimatd for all radoms effects in pois
    v=diag(c(rand_var_pois,rand_var_zi),2) ##assuming no covariance 
    
    untrans_LBS_zi_sum=sum(sols_zi[row,]*ibc)+sol_inter_zi[row,]#sum of all chromosomal effects x inbreeding coefficient 
    untrans_LBS_pois_sum=sum(sols_pois[row,]*ibc)+sol_inter_pois[row,]#sum of all chromosomal effects x inbreeding coefficient 

    mu=(c(untrans_LBS_pois_sum,untrans_LBS_zi_sum))
    
    LBS_trans=normal.zipooisson(mu,v)
    LBS_sum[[row]]=LBS_trans
  }  
  
  LBS_sum_unlist=do.call(rbind.data.frame,LBS_sum)
  
  mean<-apply(LBS_sum_unlist,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
  CI_upper<-apply(LBS_sum_unlist,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
  CI_lower<-apply(LBS_sum_unlist,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
  
  LBS_trans_mat=matrix(c(mean,CI_upper,CI_lower, paste0(ibc)),nrow = 1, ncol = 4)
  LBS_transformed[[c]]=LBS_trans_mat
  
  print(paste0("Finished ibc = ",ibc))
  
}

beep(sound=8)

LBSsum_df=do.call(rbind.data.frame,LBS_transformed)
colnames(LBSsum_df) <- c("LBS_mean","CI_upr","CI_lwr","ibc")
LBSsum_df <- sapply(LBSsum_df, as.numeric)%>%as.data.frame()

LBS_pred_sum=ggplot(data=LBSsum_df,aes(x=ibc,y=LBS_mean))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding Coefficient", y="LBS", tag = "A")+
  geom_ribbon(aes(ymin = CI_lwr, ymax = CI_upr), alpha = 0.3)+
  theme(text = element_text(size = 18), plot.title = element_text(size=13))+
  geom_point(data = male_LBS_df, (aes(x=FROH_sum_div, y=LBS)), inherit.aes = F, alpha=0.15)


LBS_pred_sum

# 
# > LBSsum_df
# LBS_mean      CI_upr       CI_lwr  ibc
# 1 11.396416602 19.41524445 6.168257e+00 0.00
# 2  3.135506335  5.71007273 1.476622e+00 0.05
# 3  0.835112367  1.81971495 2.923561e-01 0.10
# 4  0.217605382  0.59794715 4.318594e-02 0.15
# 5  0.056216923  0.19447481 6.010458e-03 0.20
# 6  0.014624130  0.06498364 7.798362e-04 0.25
# 7  0.003894285  0.02084847 8.708947e-05 0.30




# male_LBS_df%>%
#   filter(FROH_sum_div==0)%>%
  mean(male_LBS_df$LBS)


library(patchwork)

LBS_all3=(plot_spacer()+LBS_pred_sum+plot_spacer())/(plot_spacer()+forest|LBS_pred_independent+plot_spacer())

LBS_all3

# ggsave(LBS_all3,
#        file = "PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/male_LBS_all3_with_points.png",
#        width = 16,
#        height = 6)

ggsave(LBS_all3,
       file = "PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/male_LBS_all3_with_points_reorder2.png",
       width = 18,
       height = 12)



#save(LBSsum_df,LBS_male_all, file="PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/plot_dfs_male_LBS.RData")


## effect size Vs chr size 



## check if effect size is correlted with chr size 
deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  filter(!CHR %in% c("All","All_auto","X","unplaced"))

deermap$CHR=as.integer(deermap$CHR)


effect_v_size=FROH_sols_zi%>%inner_join(deermap)%>%
  select(CHR, length_Mb,solution)


zi=ggplot(effect_v_size, aes(x=length_Mb, y=solution, label=CHR)) +
  geom_point() +
  geom_smooth(method=lm, color="mediumseagreen") +
  geom_text(nudge_x = 5)+
  theme_classic()+
  theme(text = element_text(size = 18), axis.title.x = element_text(size=12))+
  labs(x="Chromosome length (Mb)", y="Slope estimate for zero inflation")+
  stat_cor(method = "pearson", label.y=1)



cor(effect_v_size$length_Mb, effect_v_size$solution)



effect_v_sizepi=FROH_sols_pois%>%inner_join(deermap)%>%
  select(CHR, length_Mb,solution)


pi=ggplot(effect_v_sizepi, aes(x=length_Mb, y=solution, label=CHR)) +
  geom_point() +
  geom_smooth(method=lm, color="chocolate3") +
  geom_text(nudge_x = 5)+
  theme_classic()+
  theme(text = element_text(size = 18), axis.title.x = element_text(size=12))+
  labs(x="Chromosome length (Mb)", y="Slope estimate for poisson")+
  stat_cor(method = "pearson")



chr_v_eff=zi+pi
chr_v_eff


top=LBS_pred_sum+forest+plot_layout(widths = c(4,5))


bottom=LBS_pred_independent+chr_v_eff+plot_layout(widths = c(4,5))

LBS_all4=top/bottom

LBS_all4

ggsave(LBS_all4,
       file = "PhD_4th_yr/Inbreeding_depression_models/LBS/Chapter_plots/male_LBS_all4.png",
       width = 14,
       height = 12)



