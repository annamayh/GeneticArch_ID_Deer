library(RODBC)
library(tidyverse)



db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.1.accdb" #open connection
con<-odbcConnectAccess2007(db)

female_LBS<-sqlFetch(con, "sys_LBS+LRS") #built in query - Calculate LBS and LRS for dead hinds 
Dlife<-sqlFetch(con, "sys_DLife") #built in query age in deer year

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthYear, MumCode,DeathType)

odbcClose(con)


setwd("H:/")

FROH<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,CHR,KB) %>% dplyr::rename(Code=IID)%>%filter(nchar(Code)==5)

KB_perLG<-FROH%>%dplyr::group_by(Code, CHR)%>%
  dplyr::summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 

deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  dplyr::filter(!CHR %in% c("All","All_auto","X","unplaced"))

froh_per_chr<-plyr::join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  dplyr::select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)

colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")



## LBS for females older than 14 but that are still alive 
# have most likely finished reproducing now
female_still_alive=Dlife%>%
  filter(DeerYear==2021)%>% ##most up to date year
  filter(Age>=14)%>%
  select(Code, Age)%>%
  inner_join(female_LBS)%>%
  inner_join(froh_per_chr)%>%
  left_join(life)%>%
  filter(if_any(DeathType, is.na)) %>%##10 still alove
  select(-DeathType, -Age)
## checked as of 07/2023 that all 10 females are still alive!


female_LBS_df=female_LBS%>%
  inner_join(froh_per_chr)%>%
  left_join(life)%>%
  filter(DeathType!= "S"| is.na(DeathType))%>% #keep ids that died a natural death also filters out ids with no recorded death type 
  filter(DeathType!= "A" | is.na(DeathType)) %>% ##sample size reduced from 1498 to 1070
  filter(DeathType!= "D"| is.na(DeathType) )%>%
  filter(BirthYear<=2006)%>% ## remove females younger than 14 as may not have reached full LBS yet (also removes dead ids in this time to avoid biasing data)
  select(-DeathType)%>%
  rbind(female_still_alive)%>%
  na.omit()


write.table(female_LBS_df,
            file = "PhD_4th_yr/Inbreeding_depression_models/LBS/female_LBS_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 


