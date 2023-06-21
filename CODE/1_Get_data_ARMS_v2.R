#############################################################################################################
##   Compile ARMS data from 4 datasets and create dataframes for analysis.
##   1.Pelagic genetic. 2018-2019. Mendexa. A. Lanzen.
##   2.Pelagic traditional. 2018-2019. Mendexa. I. Muxika.
##   3.Benthic genetic. 2013-2014. Lekeitio, Zumaia, Pasaia. Pearman et al. 2020 (A. Lanzen).
##   4.Benthic traditional. 2013-2014. Lekeitio, Zumaia, Pasaia. Romain et al. 2019 (I. Muxika).
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es 
##   Last updated: 19/06/2023
##############################################################################################################

#load libraries
library (tidyverse) 
library (readxl) 
library (vegan)

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 1. Genetic pelagic and benthic ARMS----   
# 1.0 read data
df <- read.delim(".././DATA/GENETIC/SWARM_table_final.tsv", row.names = 1)#556 otus

# 1.0.1 transpose
df_t <- as.data.frame(t(df[-39])) #remove classification col

# 1.0.2 calculate relative abundance of ASVs
summary (df_t)
df_ra <- decostand (df_t, method="total")
df_ra_t <- as.data.frame(t(df_ra))
df_ra_t$SWARM <- row.names(df_ra_t)
df_ra_t$classification <- df$classification

#check
table (row.names (df_ra_t) == row.names(df))
names (df_ra_t)

# 1.0.3 remove 1000-2000 size fraction only present in Mendexa
df_ra_t <-  select(df_ra_t, -Pelagic01_1000uM)


# 1.1 add metadata to ASV dataframe
samples <- read.csv(".././DATA/GENETIC/samples.csv", sep=";")#556 otus

#1.2 manipulate data, wide to long and add domain, trait, site, and size variables (there is no consistent pattern in variables so we do it manually)
df1 <- gather(df_ra_t, Sample, Abundance, Pelagic01_2000uM:ZUM.03_500uM)#
df1$Domain <- rep (c("Pelagic","Benthic"), times=c(2224,18348))
df1$Trait<-"Mobile"
df1[grep("Placa",df1$Sample),"Trait"]<-"Sessile"
df1$Site <- rep (c("Mendexa-pelagic","Lekeitio", "Pasaia", "Zumaia"), times=c(2224,6116,6116,6116))
df1$Size <- rep (c(">2000","500-1000","90-500", "Sessile",
                   "Sessile","90-500","500-1000","Sessile","90-500","500-1000",">2000","Sessile","90-500","500-1000",">2000",
                   "Sessile","90-500","500-1000","Sessile","90-500","500-1000",">2000","Sessile","90-500","500-1000",">2000",
                   "Sessile","90-500","500-1000",">2000", "Sessile","90-500","500-1000",">2000", "Sessile", "90-500","500-1000"), 
                 times=c(556,556,556,556,
                         556,556,556,556,556,556,556,556,556,556,556,
                         556,556,556,556,556,556,556,556,556,556,556,
                         556,556,556,556,556,556,556,556,556,556,556))  

# 1.2.1 add Phylum and Class columns
df1 <- df1  %>% rowwise() %>% mutate (Phylum = unlist(strsplit(classification, split = ";"))[5], Class = unlist(strsplit(classification, split = ";"))[6])

# 1.2.2 add Replicates column
df1$Repl <- rep (c("1","2","3","1","2","3","1","2","3"), times=c(3892,2224,2224,1668,2224,2224,2224,2224,1668))

#1.2.3 add lat and long
df1$Lat <- rep(c(43.5,43.5,43.4,43.4), times=c(2224,6116,6116,6116))
df1$Long <- rep(c(-2.45,-2.50,-1.93,-2.23), times=c(2224,6116,6116,6116))

# 1.2.4 do replicate means by site    
df1 <- df1   %>% group_by(Site,Lat,Long, Domain, Trait, Size, SWARM, Phylum, Class) %>% summarize(Abundance = mean(Abundance, na.rm = TRUE)) 

# 1.2.5 factor Domain
df1$Domain<- factor(df1$Domain, levels = c("Pelagic","Benthic"))
df1$Method <- "Metabarcoding" #add method column

# 1.2.6 rename
df_gen <- data.frame (df1)

# 2. Photo-analysis pelagic and benthic ARMS----
df <- data.frame(read_excel(".././DATA/PHOTO_ANALYSIS/Pelagic_Benthic/Resultados CPCe  2013-2019 analisis v1.xlsx", sheet = "Aukeratuak1", skip = 1))

# 2.1 manipulate data 
colnames(df)[1:6] <- c("Phylum","Class","Order", "Family","Genus", "Species") # add column names
df <- df[-c(1,34,35,36),] # remove algae and other useless rows
df1 <- df %>% filter (Phylum=="Porifera") %>% summarise_at(vars(7:66), sum, na.rm = TRUE) # do porifera sum, Muxika pers comm.
df1 <- cbind (df [28,1:6], df1) #join to main dataframe
df2 <- rbind (df[df$Phylum!="Porifera",],df1)
df2[30, 2] <- "Porifera" # rename
df2[30, 6] <- "Porifera" # rename
df2[29, 2] <- "Unidentified" # rename
df2 <- df2 [c(1:28,30,29),] #reorder
df2[,7:66]  <- as.data.frame(lapply(df2[,7:66], function(x) x / sum(x)))*100 # do percentages 
colSums(df2[,7:66]) # check dataframe
df2 <- df2 %>% select (-c(Lekeitio_2_1a, Pasaia_1_1a,Pasaia_2_1a,Pasaia_3_1a, Zumaia_1_1a, Zumaia_2_1a)) # remove cols with all values NA
colSums(df2[,7:60]) #check dataframe

# 2.2 wide to long and add Domain, Trait, Site and Face manually (not consistent pattern)
df3 <- gather(df2, Sample, Cover, ARMSA_P1_B:Zumaia_3_8b)
df3$Domain <- rep (c("Pelagic","Benthic"), times=c(180,1440))
df3$Trait<-"Sessile"
df3$Site <- rep (c("Mendexa-pelagic","Lekeitio", "Pasaia", "Zumaia"), times=c(180,510,450,480)) 
df3$Face <- rep (c("B","T","B","T","B","T",#men
                   "T", "B","T","B","T","B","B","T","B","T","B","T","B","T","B","T","B",#lek
                   "B","T", "B","T","B","B","T","B","T","B", "B","T","B","T","B",#pas
                   "B","T", "B","T","B","B","T","B","T","B", "T","B","T","B","T","B"), each=30) 

df3$Plate <- rep (c("1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8",
                    "1","4", "8"), times=c(60,60,60,
                                           60,60,60,
                                           30,60,60,
                                           60,60,60,
                                           30,60,60,
                                           30,60,60,
                                           30,60,60,
                                           30,60,60,
                                           30,60,60,
                                           60,60,60)) 

# 2.3 add lat and long
df3$Lat <- rep(c(43.5,43.5,43.4,43.4), times=c(180,510,450,480))
df3$Long <- rep(c(-2.45,-2.50,-1.93,-2.23), times=c(180,510,450,480))




# 2.4 replicate means by site and clean data    
df4 <- as.data.frame(df3 %>% group_by(Site, Lat, Long, Domain, Trait, Plate, Face,Species, Phylum, Class) %>% summarize(Cover = mean(Cover, na.rm = TRUE)))
df4$Domain<- factor(df4$Domain, levels = c("Pelagic","Benthic")) #factor Domain
df4$Cover <- as.numeric (df4$Cover) #convert to numeric
df4$Method <- "Photo-analysis" #add method
df4$Domain<- factor(df4$Domain, levels = c("Pelagic","Benthic")) #level Domain

# 2.5 rename photo-analysis dataframe
df_photo <- df4

# save photo and genetic dataframes
save(df_gen,df_photo, file = "./dat_arms.RData")





