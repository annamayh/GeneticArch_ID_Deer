library(tidyverse)
library(RODBC)




db<-"C:\\Users\\s1881212\\Documents\\Deer_database_2022/RedDeer2.05.1.accdb" #open connection
con<-odbcConnectAccess2007(db)

mum_birthyr<-sqlFetch(con, "tbllife")%>%select(Code,BirthYear)%>%rename(mum_birthyear=BirthYear,MumCode=Code)
#mum_stat<-sqlFetch(con, "sys_HindStatusAtConception") 
birth_wt<-sqlFetch(con, "sys_BirthWt")%>%dplyr::select(BirthWt,Code)

life<-sqlFetch(con, "tbllife") %>% 
  dplyr::select(Code, BirthDay,BirthMonth, BirthYear,DeathDay,DeathMonth,DeathYear, Sex, MumCode, DeathType)
ped<-sqlFetch(con, "sys_Pedigree")%>%dplyr::select(Code, MumCode,Sire)

cen=sqlFetch(con, "sys_LastCensusSighting")

odbcClose(con)

setwd("H:/")
mum_stat=read.csv("PhD_4th_yr/Inbreeding_depression_models/hind_status_at_conception_download072022.csv", header = T, stringsAsFactors = F)%>%
 select(-CalfBirthYear)%>%rename(MumCode=Mum, Code=Calf)

day_seq=read.csv("PhD_4th_yr/Inbreeding_depression_models/date_seq.csv", header = TRUE)
head(day_seq)
DOBs=life%>%select(Code, BirthDay, BirthMonth)%>%
  left_join(day_seq)%>% ## sequences of day starting at 30th April 
  select(Code, Day_seq)


##################################
#### calculating juvenile survival###
######################################

##based on last census sighting 
year=life

census_sep=cen%>%separate(LastSeen,into=c("YearLast","MonthLast","DayLast"),sep="-")%>%
  dplyr::select(Code, YearLast)

census_sep$YearLast=as.numeric(census_sep$YearLast)
year_plus_census=year%>%right_join(census_sep)%>%
  mutate(yrs_since_last_Seen=(2023-YearLast))

year_plus_census_check=year_plus_census%>%filter(is.na(DeathYear))%>%filter(!is.na(BirthYear))%>%
  filter(BirthYear==YearLast)

## for all ids with missing death years fill in estimated death year when havent been seen for at least 3 yrs 
year_plus_census_filt=year_plus_census%>%filter(is.na(DeathYear))%>%
  mutate(DeathYear_est=YearLast+3)%>%
  filter(yrs_since_last_Seen>1)%>% #has to be more than 3 yrs scince its been seen to have an estimated death
  filter(DeathYear_est<2023)%>% 
  dplyr::select(Code,DeathYear_est,YearLast)

year_plus_est_death=year%>%left_join(year_plus_census_filt) ##adding in an estimated death yr and year last seen for ids w/ no death year (not including ids born in the last 3 yrs )


year_plus_est_death$DeathMonth=as.numeric(year_plus_est_death$DeathMonth)
year_plus_est_death$DeathYear=as.numeric(year_plus_est_death$DeathYear)
year_plus_est_death$DeathYear_est=as.numeric(year_plus_est_death$DeathYear_est)
year_plus_est_death$BirthYear=as.numeric(year_plus_est_death$BirthYear)
year_plus_est_death$YearLast=as.numeric(year_plus_est_death$YearLast)

# year_plus_est_death_check=year_plus_est_death%>%filter(is.na(DeathYear))%>%filter(!is.na(BirthYear))%>%
#   filter(BirthYear==YearLast)

year_juve=year_plus_est_death%>%
  mutate(juvenile_survival = case_when(
       ((is.na(DeathYear)) & (BirthYear==YearLast))~"0",

    ## when there is no recorded death year use the estimated death year (which is when id has not been seen for 3 years...)
    # not sure any juveniles will actually be recorded as 0 in this case as theyre killed off after 3 year anyway 
      (is.na(DeathYear) & (DeathYear_est<(BirthYear+2)))~"0",
      #if died in same year born then juvenile surv=0
      (DeathYear==BirthYear) ~ "0", 
      #if died the year after it was born juvenle surv=0
      (DeathYear==(BirthYear+1)) ~ "0", 
      #if died before May in the second year of life then did not survival to juvenile (cut off is May)
     (DeathYear==(BirthYear+2))&(DeathMonth<5) ~ "0", 
      #any ids with no recorded death month and a death year 2 yrs after birth is given NA as we dont know if id died before or after 2nd birhday (only 5 ids)
     (DeathYear==(BirthYear+2))&(is.na(DeathMonth)) ~ "NA", 
      
      ) )%>%
  # all ids that didnt get a 0 for juvenile survival gets a 1 for successfully surviving 
  mutate(juvenile_survival =replace_na(juvenile_survival,"1"))
  

table(year_juve$DeathType)
 
juvenile_surv=year_juve%>%filter(DeathType!= "S"| is.na(DeathType))%>% #filter out ids that were shot
  filter(DeathType!= "A" | is.na(DeathType)) %>% #%>% ##filter out ids that were an accidental death
  filter(DeathType!= "D" | is.na(DeathType))%>% ## keeping NA death type
  dplyr::select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv$juvenile_survival)

# want to keep ids that were shot but also made it to adulthood 
juvenile_surv_shot=year_juve%>%filter(DeathType== "S"&juvenile_survival=="1")%>%
  select(Code, BirthYear, Sex, MumCode, juvenile_survival)

table(juvenile_surv_shot$juvenile_survival)

## adding ids that survived to age 2 then were shot to the df
juvenile_surv=juvenile_surv%>%rbind(juvenile_surv_shot)
table(juvenile_surv$juvenile_survival)


# ###########################
# ## read in FROH values ####
############################

FROH<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom", header=T, stringsAsFactors = F)%>%
  select(IID,CHR,KB) %>% rename(Code=IID)%>%filter(nchar(Code)==5)

KB_perLG<-FROH%>%group_by(Code, CHR)%>%
  summarise(KB_chr=sum(KB))%>%
  ungroup %>% 
  complete(Code, CHR, fill = list(KB_chr = 0)) #completes for all chromosomes 

deermap <- read.csv("PhD_3rdYR/Data_files/Genome_assembly_mCerEla1.1.csv", header = T, stringsAsFactors = F)%>%
  filter(!CHR %in% c("All","All_auto","X","unplaced"))

froh_per_chr<-plyr::join(KB_perLG,deermap)%>% mutate(chr_froh=KB_chr/length_Kb)%>%
  select(-length_Mb,-length)%>% reshape2::dcast(Code~CHR) %>% 
  mutate(FROHsum = rowSums(.[2:34]))%>%mutate(FROH_sum_div=FROHsum/33)

colnames(froh_per_chr) <- c("Code", paste0("FROH_chr", 1:33),"FROHsum","FROH_sum_div")

FROH_full<-read.table("PhD_4th_yr/2023_ROH_search/2021_sleuthed_052023.hom.indiv", header=T, stringsAsFactors = F)%>%
  dplyr::select(IID,KB) %>% dplyr::rename(Code=IID)%>%mutate(FROH=KB/2591865)%>%filter(nchar(Code)==5)

#combine all datasets to make a dataframe 

juvenile_surv_df=juvenile_surv%>%
  left_join(mum_stat)%>%
  left_join(birth_wt)%>%
  left_join(mum_birthyr)%>%
  mutate(mum_age=BirthYear-mum_birthyear)%>%
  inner_join(froh_per_chr)%>%## only joining when ids have FROH values and juvenile surv values
  mutate(mum_age_sq=mum_age^2)%>%
  select(-mum_birthyear)%>%#dont need mum birth yr in df anymore
  filter(!Sex=="3")%>%
  left_join(DOBs)
  
##some NAs in mum stat and birth weight


write.table(juvenile_surv_df,
            file = "PhD_4th_yr/Inbreeding_depression_models/survival/AA_juvenile_survival_df.txt",
            row.names = F, quote = F, sep = ",",na = "NA") #saving tables as txt file 

