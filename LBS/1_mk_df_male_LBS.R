## MALE LBS ###
library(RODBC)
library(tidyverse)


# Set up dataframe #

db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.1.accdb" #open connection
con<-odbcConnectAccess2007(db)
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, BirthYear,MumCode,Sire)
life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code,DeathYear,Sex,DeathType,BirthYear)
odbcClose(con)

sire_count=ped%>%dplyr::select(Code,Sire)%>%
  group_by(Sire)%>%tally()%>% ##grouping by sire and counting up ids sired by id
  rename(Code=Sire, LBS=n)

dead_males=life%>%
  filter(Sex=="2")%>% #select males
  filter(DeathType!= "S")%>% #filter out ids that were shot so we only keep ids that died a natural death 
  filter(DeathType!= "A" ) %>% ##
  filter(DeathType!= "D" )

all_male_LBS=left_join(dead_males,sire_count)%>% #before filter 1251, after956
  filter(BirthYear<=2008)
all_male_LBS[is.na(all_male_LBS)] <- 0 #males with no recorded sired ids given 0 LBS


# Read in ROH per chr #

setwd("H:/")
###need to make this with updated map positions
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

# Combine dataframes #

male_LBS_df=all_male_LBS%>%
  left_join(froh_per_chr)%>%left_join(ped)%>%
  dplyr::select(-Sire, -DeathYear, -Sex, -DeathType)%>%
  na.omit()


write.table(male_LBS_df,
            file = "PhD_4th_yr/Inbreeding_depression_models/LBS/Male_LBS_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 

