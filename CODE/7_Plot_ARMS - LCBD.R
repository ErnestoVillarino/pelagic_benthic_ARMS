#############################################################################################################
##   Local Contribution to Betadiversity in pelagic and benthic habitats
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################

# load libraries
library (adespatial) #LCBD
library (tidyverse) #manipulate data
library (cowplot) #arrange plots

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")
df_gen$Site  <- factor(df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia")) 
df_photo$Site<- factor(df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1. Manipulate genetic and photo-analysis datasets long to wide format and do means according to group
df1  <-  df_gen    %>% 
  group_by(Site,SWARM,Size, Domain) %>% 
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  spread(SWARM,Abundance)

df2  <-  df_photo   %>% 
  group_by(Site,Species,Face, Plate, Domain) %>% 
  summarize(Cover = mean(Cover, na.rm = TRUE)) %>% 
  spread(Species,Cover) %>% 
  select (-Unidentified) 

# 1. Local contribution to beta-diversity-----

# 1.1 LCBD gen
res1 <-  beta.div(df1[,4:559], "hellinger", nperm=999)
df1$LCBD <- as.vector(res1[["LCBD"]])
df1$Method <- "Metabarcoding"

# 2. LCBD photo-analysis
res2 <-  beta.div(df2[,5:33], "hellinger", nperm=999)
df2$LCBD <- as.vector(res2[["LCBD"]])
df2$Method <- "Photo-analysis"

# plot LCBD gen
summary (df1$LCBD) # take  mean

df1$var <- "Habitat"
p1 <- ggplot(df1, aes(x=Domain, y=LCBD, fill=Domain))+
  xlab("")+
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") +
  theme_bw() + facet_grid(~var+Method)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999")) +
  theme(axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text(size=16),
        axis.text.y = element_text (size=16),
        strip.text = element_text(size=16), legend.position="none",
        panel.grid.major = element_blank())  + geom_hline(yintercept=0.05882, linetype="dashed", color = "black") 

# plot LCBD photo-analysis
summary (df2$LCBD) # take  mean

df2$var <- "Habitat"
p2 <- ggplot(df2, aes(x=Domain, y=LCBD, fill=Domain))+
  xlab("")+
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") +
  theme_bw() + facet_grid(~var+Method)+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999")) +
  theme(axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text(size=16),
        axis.text.y = element_text (size=16),
        strip.text = element_text(size=16), legend.position="none",
        panel.grid.major = element_blank())  + geom_hline(yintercept=0.05882, linetype="dashed", color = "black") 

#arrange plots
plot_grid(p1,p2, labels = c('a', 'b'), label_size = 22, ncol=2, nrow = 2)

#save
ggsave(".././FIGURES/F8ab.png",  width=14, height=14, dpi=300, bg = "white")

