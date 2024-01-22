
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("H:/")
load(file="PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_model_output24072023.RData")

birth_wt_df=read.table("PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_df.txt", sep=",", header=T)
birth_wt_df_na_rm=birth_wt_df%>%
  na.omit()

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



forest=ggplot(data=FROH_sols, aes(x=CHR, y=solution, ymin=CI_lower, ymax=CI_upper)) +
  geom_pointrange(colour="grey50") + #plots lines based on Y and lower and upper CI
  geom_hline(yintercept=0, lty=2) +  # black line is 0
  coord_flip() +  # flip coordinates (puts labels on y axis)
  labs(x="Chromosome", y="Posterior estimate + CI", title = "Variation in chromosome-specific inbreeding effects") +
  theme_classic()+  # use a white background                 
  theme(legend.position = "none")+
  ## line of FROHsum estimate
  annotate("segment", x = 0, xend = 34, y = FROH_sum_sol, yend = FROH_sum_sol, colour = "red", alpha=0.6)+
    ##add in estimate of FROHsum
  annotate("pointrange", x = 34, y = FROH_sum_sol, ymin = FROH_sum_hu_upr, ymax = FROH_sum_hu_lwr,
           colour = "red", linewidth = 1, alpha=0.5, size=0.2)+
  ##significance of FROHsum
  geom_text(aes(x=34.5, y=FROH_sum_sol, label="***"), colour="red", alpha=0.5)+
  expand_limits(x = 35)+
  theme(text = element_text(size = 18), plot.title = element_text(size=13), axis.text.y = element_text(size=10))



forest

###############################################################################################################
#########################################################################################
##########################################################################################################


# model Sex + AgeHrs + MotherStatus + mum_age+mum_age_sq+FROHsum+(mm chr FROH)+mumcode(random)+birthyr(random)

#read in table to see what averages are
# birth_wt_df=read.table("PhD_4th_yr/Inbreeding_depression_models/birth_weight/birth_wt_df2023.txt", sep=",", header=T)
# birth_wt_df_na_rm=birth_wt_df%>%
#   dplyr::select(-mum_birthyear)%>%na.omit()

##check mean birth weight
mean(birth_wt_df_na_rm$CaptureWt) ##this is what we are predicting 
#7.097226
#mean(birth_wt_df_na_rm$BirthWt)
#6.483708


## finding means coefficients from df used 
Sex=0.5#1=female 2=male intercept is females 
AgeHrs=median(birth_wt_df_na_rm$AgeHrs)#median is 24 hrs 
mum_age=mean(birth_wt_df_na_rm$mum_age) #mean mum age in dataset
mum_age_sq=mum_age^2
day_seq=33


##also choose quantiles 
ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)
 

pred_ibcs=list()

for(v in 1:length(ibc_qua)){

    ibc=ibc_qua[v] #picking inbreeding coeff from quantiles vector
   # FROHsum=(ibc*33)##choose inbreeding coefficient and * by the number of chr (which is what FROHsum is)
    #FROHchrs=ibc
    
    bw_pred_list=list()
    
    for (i in 1:33){
      bw_pred=summary_table[1,1]+#intercept
        (Sex *(summary_table[2,1]))+ 
        (AgeHrs*(summary_table[3,1])) + #Assuming calf is 0 hours old 
        (mum_age*(summary_table[8,1]))+
        (mum_age_sq*(summary_table[9,1]))+
        (day_seq*summary_table[10,1])+
        (ibc*(as.numeric(FROH_sols[i,1])))
      bw_pred_list[i]=bw_pred
    }
    
    bw_pred_chr=do.call(rbind.data.frame,bw_pred_list)%>%
      tibble::rownames_to_column("CHR")%>%
      add_column(paste0(round(ibc, digits = 4)))
     
    colnames(bw_pred_chr) <- c("CHR", "bw_prediction","ibc")
    
    pred_ibcs[[v]]=bw_pred_chr
    
}




pred_ibcs_all=do.call(rbind.data.frame,pred_ibcs)

pred_ibcs_all$CHR=as.numeric(pred_ibcs_all$CHR)


chrind=pred_ibcs_all%>%arrange(CHR)%>%
  ggplot(aes(x=as.numeric(ibc),y=bw_prediction, group=as.factor(CHR), colour=as.factor(CHR)))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x=" Inbreeding coefficient", y="Predicted birth weight (kg)", title="Assuming chromosome independence",colour="Chromosome" )+
    theme(text = element_text(size = 18), 
          plot.title = element_text(size=13), 
          legend.title = element_text(size=10),
          legend.text = element_text(size=9))


chrind





### now assumung indiv is inbred by same amount on all chromosomes
# plus CIs

Sex=0.5#1=female 2=male intercept is females 
AgeHrs=median(birth_wt_df_na_rm$AgeHrs)#median is 24 hrs 
mum_age=mean(birth_wt_df_na_rm$mum_age) #mean mum age in dataset
mum_age_sq=mum_age^2
day_seq=33


inter=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(matches("Interc"))

FROHsumest=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(matches("FROHs"))

sex_ests=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(matches("Sex"))
sex_ests=sex_ests*Sex

age_ests=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(("AgeHrs"))
age_ests=age_ests*AgeHrs

mumage_ests=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(("mum_age"))
mumage_ests=mumage_ests*mum_age

mumagesq_ests=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(("mum_age_sq"))
mumagesq_ests=mumagesq_ests*(mum_age^2)

day_seq_ests=as.data.frame(birth_wt_model$Sol)%>%dplyr::select(("Day_seq"))
day_seq_ests=day_seq_ests*day_seq


all_effects=as.matrix(rowSums(inter%>%cbind(sex_ests)%>%cbind(age_ests)%>%cbind(mumage_ests)%>%cbind(mumagesq_ests)%>%cbind(day_seq_ests)))


ibc_qua=c(0,0.05,0.1,0.15,0.2,0.25,0.3)


pred_nonind=list()
for(v in 1:length(ibc_qua)){
  
    ibc=ibc_qua[v]
    pred_per_it=list()
      
    
    for (row in 1:nrow(sols_full)){
    iter_chr=sum(sols_full[row,]*ibc)#sum of all chromosomal effects x inbreeding coefficient 
    iter_fixed=all_effects[row,]
    pred_per_it[[row]]=iter_chr+iter_fixed
    
      
    }
    
    bw_pred_unlist=do.call(rbind.data.frame,pred_per_it)
    
    mean<-apply(bw_pred_unlist,2,mean)#gets mean of all solutions i.e. the effect size of random effects 
    CI_upper<-apply(bw_pred_unlist,2,quantile,probs = c(0.975)) #gets upper confidence interval for all solutions 
    CI_lower<-apply(bw_pred_unlist,2,quantile,probs = c(0.025)) #gets lower CI for all solutions
    
    bw_mat=matrix(c(mean,CI_upper,CI_lower, paste0(ibc)),nrow = 1, ncol = 4)
    pred_nonind[[v]]=bw_mat
    

    }


pred_ibcs_mean=do.call(rbind.data.frame,pred_nonind) 

names(pred_ibcs_mean)<-c("Mean","upperCI","lowerCI","ibc")
pred_ibcs_mean <- sapply(pred_ibcs_mean, as.numeric)%>%as.data.frame()


birth_wt_df_na_rm_plot=birth_wt_df_na_rm%>%filter(AgeHrs<48)

chrnonind=ggplot(data=pred_ibcs_mean,aes(x=ibc,y=Mean))+
  geom_line(linewidth=1)+
  theme_bw()+
  labs(x="Inbreeding coefficient", y="Predicted birth weight (kg)", title="Assuming equal inbreeding on all chromosomes")+
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = 0.1)+
  theme(text = element_text(size = 18), 
        plot.title = element_text(size=13))+
  geom_point(data = birth_wt_df_na_rm_plot, (aes(x=FROH_sum_div, y=BirthWt)), inherit.aes = F, alpha=0.15)+
  xlim(0,0.3)

chrnonind

library(patchwork)

Bw_3in1=chrnonind+forest+chrind+plot_annotation(tag_levels = 'A')

Bw_3in1

# ggsave(Bw_3in1,
#        file = "PhD_4th_yr/Inbreeding_depression_models/birth_weight/all3_birthweight_points.png",
#        width = 17,
#        height = 6)
# 

### effect size Vs chr size


## check if effect size is correlted with chr size 
deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  filter(!CHR %in% c("All","All_auto","X","unplaced"))


effect_v_size=FROH_sols%>%inner_join(deermap)%>%
  select(CHR, length_Mb,solution)

effect_v_size$CHR=as.factor(effect_v_size$CHR)

effect_chr=ggplot(effect_v_size, aes(x=length_Mb, y=solution, label=CHR)) +
  geom_point() +
  geom_smooth(method=lm, color="red") +
  geom_text(nudge_x = 3)+
  theme_classic()+
  theme(text = element_text(size = 18))+
  labs(x="Chromosome length (Mb)", y="Slope estimate")+
  stat_cor(method = "pearson",label.x = 110)

effect_chr

cor(effect_v_size$length_Mb, effect_v_size$solution)



Bw_4in1=(chrnonind+forest)/(chrind+effect_chr)+plot_annotation(tag_levels = 'A')

Bw_4in1


ggsave(Bw_4in1,
       file = "PhD_4th_yr/Inbreeding_depression_models/birth_weight/all4_birthweight.png",
       width = 13,
       height = 12)

