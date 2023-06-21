#############################################################################################################
##   Plot pair site Bray-Curtis dissimilarity matrices.
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################

# load libraries
library (vegan)
library (tidyverse)
library (ggcorrplot)
library (patchwork) 

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")
df_gen$Site  <- factor(df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))
df_photo$Site<- factor(df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1. Dissimilarity Bray-Curtis (genetic)-----

# 1.1 Group by Site
df_gen_site <-  df_gen  %>% group_by(Site,SWARM) %>% 
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  spread(SWARM,Abundance)

# 1.1.1 Calculate Bray-Curtis dissimilarity
df_bc_gen_site <- vegdist (df_gen_site[2:557], "bray", na.rm = T)  # Bray-Curtis.
df1 <- as.matrix (df_bc_gen_site)  #convert to matrix
colnames(df1) <-c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia")
rownames(df1) <-c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia") 

# 1.2 Group by Site and Size
df_gen_size_site <-  df_gen   %>% 
  group_by(Site,SWARM,Size) %>% 
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  spread(SWARM,Abundance)

# 1.2.1 Calculate Bray-Curtis dissimilarity
df_bc_gen_size_site <- vegdist (df_gen_size_site[3:558], "bray", na.rm = T)  # Bray-Curtis.
df2 <- as.matrix (df_bc_gen_size_site) #convert to matrix
colnames(df2) <-c("Mendexa-pelagic >2000", "Mendexa-pelagic 500-1000", "Mendexa-pelagic 90-500", "Mendexa-pelagic Sessile",
                  "Lekeitio >2000", "Lekeitio 500-1000","Lekeitio 90-500","Lekeitio Sessile",
                  "Zumaia >2000", "Zumaia 500-1000","Zumaia 90-500","Zumaia Sessile",
                  "Pasaia >2000", "Pasaia 500-1000","Pasaia 90-500","Pasaia Sessile")

rownames(df2) <-c("Mendexa-pelagic >2000", "Mendexa-pelagic 500-1000", "Mendexa-pelagic 90-500", "Mendexa-pelagic Sessile",
                  "Lekeitio >2000", "Lekeitio 500-1000","Lekeitio 90-500","Lekeitio Sessile",
                  "Zumaia >2000", "Zumaia 500-1000","Zumaia 90-500","Zumaia Sessile",
                  "Pasaia >2000", "Pasaia 500-1000","Pasaia 90-500","Pasaia Sessile")

#plot gen site
p1 <- ggcorrplot(df1, type = "lower",lab = T, show.diag = T, lab_size = 5) + 
  scale_fill_distiller(limit = c(0,1), palette = "RdPu",guide = guide_colorbar(frame.colour = "black")) + 
  labs(fill = "Community\ndissimilarity") +
  theme(axis.text.x = element_text (size=16),
        axis.text.y = element_text (size=16),
        legend.title = element_text( size = 16), 
        legend.text = element_text(size=16)) + ggtitle("Metabarcoding")

# plot gen site size
p2 <- ggcorrplot(df2, type = "lower",lab = F, show.diag = T, lab_size = 5) + 
  scale_fill_distiller(limit = c(0,1), palette = "RdPu",guide = guide_colorbar(frame.colour = "black")) + 
  labs(fill = "Community\ndissimilarity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13),
        axis.text.y = element_text (size=13),
        legend.position = "none") + ggtitle("Metabarcoding")


# 2. Dissimilarity photo-analysis------

# 2.1 Group by site
df_photo_site <- df_photo %>% 
  group_by(Site,Species) %>% 
  summarize(Cover = mean(Cover, na.rm = TRUE)) %>% 
  spread(Species,Cover) %>% select (-Unidentified) 

# 2.1.1 Calculate Bray-Curtis dissimilarity
df_bc_photo_site <- vegdist (df_photo_site[2:30], "bray", na.rm = T)  # Bray-Curtis.
df3 <- as.matrix (df_bc_photo_site) #convert to matrix
colnames(df3) <-c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia")
rownames(df3) <-c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia")

# 2.2 Group by site, face and plate photoitional sites, face, plate
df_photo_site_face <- df_photo  %>% 
  group_by(Site,Species,Face,Plate) %>% 
  summarize(Cover = mean(Cover, na.rm = TRUE)) %>% 
  spread(Species,Cover) %>% 
  select (-Unidentified) 

# 2.2.1 Calculate Bray-Curtis dissimilarity
df_bc_photo_site_face <- vegdist (df_photo_site_face[4:32], "bray", na.rm = T)  # Bray-Curtis.
df4 <- as.matrix (df_bc_photo_site_face) #convert to matrix
colnames(df4) <-c("Mendexa-pelagic B1", "Mendexa-pelagic B4", "Mendexa-pelagic B8", "Mendexa-pelagic T1", "Mendexa-pelagic T4", "Mendexa-pelagic T8",
                  "Lekeitio B1", "Lekeitio B4","Lekeitio B8","Lekeitio T1", "Lekeitio T4","Lekeitio T8",
                  "Zumaia B1", "Zumaia B4","Zumaia B8","Zumaia T1", "Zumaia T4","Zumaia T8",
                  "Pasaia B1", "Pasaia B4","Pasaia B8","Pasaia T4", "Pasaia T8")

rownames(df4) <-c("Mendexa-pelagic B1", "Mendexa-pelagic B4", "Mendexa-pelagic B8", "Mendexa-pelagic T1", "Mendexa-pelagic T4", "Mendexa-pelagic T8",
                  "Lekeitio B1", "Lekeitio B4","Lekeitio B8","Lekeitio T1", "Lekeitio T4","Lekeitio T8",
                  "Zumaia B1", "Zumaia B4","Zumaia B8","Zumaia T1", "Zumaia T4","Zumaia T8",
                  "Pasaia B1", "Pasaia B4","Pasaia B8","Pasaia T4", "Pasaia T8")

# plot photo-analysis site
p3 <- ggcorrplot(df3, type = "lower",lab = TRUE, show.diag = T, lab_size = 5) + 
  scale_fill_distiller(limit = c(0,1), palette = "RdPu",guide = guide_colorbar(frame.colour = "black")) + 
  labs(fill = "Community\ndissimilarity") +
  theme(axis.text.x = element_text (size=16),
        axis.text.y = element_text (size=16),
        legend.title = element_text( size = 16), 
        legend.text = element_text(size=16))+ ggtitle("Photo-analysis")

# plot photo-analysis site face plate
p4 <- ggcorrplot(df4, type = "lower",lab = F, show.diag = T, lab_size = 5) + 
  scale_fill_distiller(limit = c(0,1), palette = "RdPu",guide = guide_colorbar(frame.colour = "black")) + 
  labs(fill = "Community\ndissimilarity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
        axis.text.y = element_text (size=12),
        legend.position = "none")  + ggtitle("Photo-analysis")

# arrange plots
p1+p2+p3+p4+ 
  plot_annotation(tag_levels = 'a') + 
  plot_layout (ncol=2, nrow=2, guides="collect") & 
  theme(plot.tag = element_text(size = 28))

#save
ggsave(".././FIGURES/F6.png",  width=14, height=16, dpi=300)


  

