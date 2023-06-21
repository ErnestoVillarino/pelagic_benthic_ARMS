#############################################################################################################
##   Uniqueness analysis in photo-analysis and metabarcding ARMS 
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################


## 0. load libraries
library (tidyverse)
library (ggvenn)  #venn diagram 

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")
df_gen$Site  <- factor(df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))
df_photo$Site<- factor(df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1. Unique genetic-----
# remove zeros
df1 <- df_gen[df_gen$Abundance != 0, ] #remove zeros

# filter data according to group and count unique ASVs
df1 %>% filter(Domain=="Pelagic") %>% summarize(n = n_distinct(SWARM)) #158 ASVs pelagic
df1 %>% filter(Domain=="Benthic") %>% summarize(n = n_distinct(SWARM)) #530 ASVs benthic
df1 %>% filter(Site=="Zumaia") %>% summarize(n = n_distinct(SWARM))    #349 ASVs zumaia
df1 %>% filter(Site=="Lekeitio") %>% summarize(n = n_distinct(SWARM))  #330 ASVs lekeitio
df1 %>% filter(Site=="Pasaia") %>% summarize(n = n_distinct(SWARM))    #333 ASVs pasaia
df1 %>% filter(Site=="Mendexa-pelagic") %>% summarize(n = n_distinct(SWARM)) #158 ASVs mendexa
df1 %>% filter(Size==">2000") %>% summarize(n = n_distinct(SWARM))    #326 ASVs >2000
df1 %>% filter(Size=="500-1000") %>% summarize(n = n_distinct(SWARM))    #423 ASVs 500-1000
df1 %>% filter(Size=="90-500") %>% summarize(n = n_distinct(SWARM))    #481 ASVs 90-500
df1 %>% filter(Trait=="Sessile") %>% summarize(n = n_distinct(SWARM))    #365 ASVs sessile
df1 %>% filter(Trait=="Mobile") %>% summarize(n = n_distinct(SWARM))    #540 ASVs sessile


#domain specific set of ASVs
pelagic_unique <- data.frame(SWARM=unique (df1$SWARM [df1$Domain=="Pelagic"])[!unique (df1$SWARM [df1$Domain=="Pelagic"]) %in% unique (df1$SWARM [df1$Domain=="Benthic"])]) #26
benthic_unique <- data.frame(SWARM=unique (df1$SWARM [df1$Domain=="Benthic"])[!unique (df1$SWARM [df1$Domain=="Benthic"]) %in% unique (df1$SWARM [df1$Domain=="Pelagic"])]) #398

sum (unique(df1$SWARM [df1$Domain=="Benthic"]) %in% unique (df1$SWARM [df1$Domain=="Pelagic"])) #132
sum (unique(df1$SWARM [df1$Domain=="Pelagic"]) %in% unique (df1$SWARM [df1$Domain=="Benthic"])) #132

#domain specific set of species with merge
merge_pelagic <- pelagic_unique %>% left_join (df1 %>% filter(Domain=="Pelagic"), by="SWARM")
merge_benthic <- benthic_unique %>% left_join (df1 %>% filter(Domain=="Benthic"), by="SWARM")

#check phylum
g <- merge_pelagic [! duplicated (merge_pelagic$SWARM),]
unique (merge_pelagic$Phylum)

r <- merge_benthic [! duplicated (merge_benthic$SWARM),]
unique (merge_benthic$Phylum)


#2. Unique photo-analysis------

# remove rows with 0 abundance
df2 <- df_photo[df_photo$Cover != 0, ]


# filter data according to group and count unique species
df2 %>% filter(Domain=="Pelagic") %>% summarize(n = n_distinct(Species)) #17 sp pelagic
df2 %>% filter(Domain=="Benthic") %>% summarize(n = n_distinct(Species)) #18 sp benthic
df2 %>% filter(Site=="Zumaia") %>% summarize(n = n_distinct(Species))    #10 sp zumaia
df2 %>% filter(Site=="Lekeitio") %>% summarize(n = n_distinct(Species))  #9 sp lekeitio
df2 %>% filter(Site=="Pasaia") %>% summarize(n = n_distinct(Species))    #10 sp pasaia
df2 %>% filter(Site=="Mendexa-pelagic") %>% summarize(n = n_distinct(Species)) #17 sp mendexa
df2 %>% filter(Plate=="1") %>% summarize(n = n_distinct(Species))    #20 sp p1
df2 %>% filter(Plate=="4") %>% summarize(n = n_distinct(Species))    #19 sp p4
df2 %>% filter(Plate=="8") %>% summarize(n = n_distinct(Species))    #21 sp p8
df2 %>% filter(Face=="B") %>% summarize(n = n_distinct(Species))    #23 sp top
df2 %>% filter(Face=="T") %>% summarize(n = n_distinct(Species))    #23 sp bottom



# 3. Venn diagram Two dimensions----- 

#filter data by habitat
df_pelagic <- df1 %>% filter(Domain=="Pelagic") 
df_benthic <- df1 %>% filter(Domain=="Benthic") 

#plot unique pelagic, unique benthic, shared
p1 <- ggvenn(list(Pelagic=df_pelagic$SWARM, Benthic=df_benthic$SWARM),
       stroke_size = 1,
       text_size = 7,
       set_name_size = 7,
       auto_scale = T)


#save
ggsave(".././FIGURES/F8c.png",  width=10, height=10, dpi=300, bg="white") 


  
