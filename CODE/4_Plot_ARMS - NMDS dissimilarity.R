#############################################################################################################
##   Plot ARMS NMDS on Bray-Curtis dissimilarity distances
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################

#load libraries
library (tidyverse) 
library (vegan)
library (patchwork)

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")

# 1. NMDS genetics-----

# 1.1 data manipulation, long to wide format 
df_gen <- df_gen %>% 
  select (Domain, Site, Size, Trait, SWARM, Abundance)  %>% 
  group_by (Domain, Site, Size, Trait, SWARM)  %>% 
  spread (SWARM,Abundance)

# 1.2 Bray Curtis distance on SWARMs
df_bc_gen <- vegdist (df_gen[5:560], "bray", na.rm = T)  # Bray-Curtis.

# 1.3 Calculate the "coordinates" of each sample in the NMDS space
NMDSres_gen<- metaMDS(df_bc_gen)
NMDSres_gen$stress # 0.1298309

# 1.4 Prepare data to plot (Use the scores function from vegan to extract the site scores and convert to a data.frame)
data.scores_df_gen <- as.data.frame(NMDSres_gen[["points"]]) 

# 1.5 add variables to dataframe
data.scores_df_gen$Domain <- df_gen$Domain 
data.scores_df_gen$Trait <- df_gen$Trait
data.scores_df_gen$Site <- df_gen$Site
data.scores_df_gen$Size <- df_gen$Size
data.scores_df_gen$Method <- "Metabarcoding"
data.scores_df_gen$Site<- factor(data.scores_df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1.6 plot nmds Bray-Curtis dissimilarity (genetic)
p1 <- ggplot() +  geom_point(data=data.scores_df_gen,aes(x=MDS1,y=MDS2, colour=Site, shape = Size),size=5) + 
  theme_bw() + facet_grid(~Method) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999"))+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=16),
        axis.text.y = element_text (size=16), 
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.position = "right",legend.box="vertical",
        legend.background = element_rect(colour ="black", linewidth =1, linetype="solid"),
        panel.grid.major = element_blank())  + annotate("text", x = -0.30, y = -0.21, label = "Stress= 0.12", size=6)


#2. NMDS photo-analysis------

#2.1 Manipulate data, long to wide and remove unidentified taxa
df_photo <- df_photo %>% 
  select (Site, Species, Cover, Face, Domain, Plate)  %>% 
  group_by (Site, Species, Face, Domain, Plate)  %>% 
  spread (Species,Cover) %>% 
  select (-Unidentified) 

# 2.2 rename variables
df_photo$Face [df_photo$Face=="B"]  <- "Bottom_face"
df_photo$Face [df_photo$Face=="T"]  <- "Top_face"
 
# 2.3 Bray Curtis distance on species
df_bc_photo <- vegdist (df_photo[5:33], "bray", na.rm = T)  # Bray-Curtis.

# 2.4 Calculate the "coordinates" of each sample in the NMDS space
NMDSres_photo <- metaMDS(df_bc_photo)
NMDSres_photo$stress #0.057

# 2.5 Prepare data to plot (Use the scores function from vegan to extract the site scores and convert to a data.frame)
data.scores_df_photo <- as.data.frame(NMDSres_photo[["points"]])   

# 2.6 add variables to dataframe
data.scores_df_photo$Domain <- df_photo$Domain 
data.scores_df_photo$Site <- df_photo$Site 
data.scores_df_photo$Face <- df_photo$Face 
data.scores_df_photo$Method <- "Photo-analysis"
data.scores_df_photo$Site<- factor(data.scores_df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

#2.7 plot nmds Bray - Curtis distances (photo-analysis)
p2 <- ggplot() +  geom_point(data=data.scores_df_photo,aes(x=MDS1,y=MDS2, colour=Site, shape=Face),size=5) + 
  theme_bw() + facet_grid(~Method) + guides(color = guide_legend(order=1),shape = guide_legend(order=2))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999"))+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=16),
        axis.text.y = element_text (size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.position="right",legend.box="vertical",
        legend.background = element_rect(colour ="black", size=1, linetype="solid"),
        panel.grid.major = element_blank())  + annotate("text", x = -0.46, y = -0.14, label = "Stress= 0.06", size=6)


p1/p2+ plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 28)) # NMDS

ggsave(".././FIGURES/F5.png",  width=10, height=16, dpi=300)






