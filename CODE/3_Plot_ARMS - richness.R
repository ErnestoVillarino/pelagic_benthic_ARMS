#############################################################################################################
##   Plot ARMS ASVs and species richness by site, size, plate face etc...
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################

#load libraries
library (ggplot2)
library (cowplot) 
library (tidyverse)


#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")
df_gen$Site  <- factor(df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))
df_photo$Site<- factor(df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1. ASVs richness---- 

# 1.1 calculate richness by group and do dataframe to plot (p_a)
df <- df_gen[df_gen$Abundance != 0, ] #remove zeros

df %>%
  group_by(Site) %>%
  summarize(count_distinct = n_distinct(SWARM))

df %>%
  group_by(Size) %>%
  summarize(count_distinct = n_distinct(SWARM))

df %>%
  group_by(Domain) %>%
  summarize(count_distinct = n_distinct(SWARM))

df %>%
  group_by(Trait) %>%
  summarize(count_distinct = n_distinct(SWARM))

#1.2 calculate richness by many factors
df <- df %>% group_by (Site, Size, Domain, Trait) %>% summarise(freq=n()) #
df$Method <- "Metabarcoding" #add method variable
df$Site<- factor(df$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

#plot site
df$var <- "Site"
p1 <- ggplot (df, aes(x=Site, y=freq, fill=Site))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() + #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(axis.text.x = element_text(size=16),
        legend.position="none",
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.y = element_text (size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "ASV richness") 

#plot size
df$var <- "Size"
p2 <- ggplot (df, aes(x=Size, y=freq, fill=Size))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  scale_fill_manual(values=c(">2000"="#00AFBB","500-1000"="#E7B800", "90-500"="#FC4E07", "Sessile"="#999999"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(size=16),
        legend.position="none",
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.y = element_text (size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "ASV richness") 

#plot habitat
df$var <- "Habitat"
p3 <- ggplot (df, aes(x=Domain, y=freq, fill=Domain))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  scale_fill_manual(values=c("#00AFBB","#E7B800"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(size=16),
        legend.position="none",
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.y = element_text (size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "ASV richness") 

#plot trait
df$var <- "Trait"
p4 <- ggplot (df, aes(x=Trait, y=freq, fill=Trait))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  scale_fill_manual(values=c("#00AFBB","#E7B800"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(size=16),
        legend.position="none",
        axis.title.x=element_text(size=16), # remove x-axis labels
        axis.title.y=element_text(size=16), # remove y-axis labels
        axis.text.y=element_text (size=16),
        legend.text=element_text(size=16),
        strip.text= element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "ASV richness") 



# 2. Photo-analysis richness-----

# 2.1.calculate richness by group a do dataframe to plot (p_b)
df <- df_photo[df_photo$Cover != 0, ] #remove zeros

df %>%
  group_by(Site) %>%
  summarize(count_distinct = n_distinct(Species))

df %>%
  group_by(Face) %>%
  summarize(count_distinct = n_distinct(Species))

df %>%
  group_by(Domain) %>%
  summarize(count_distinct = n_distinct(Species))

df %>%
  group_by(Trait) %>%
  summarize(count_distinct = n_distinct(Species))

#2.2 calculate richness by many factors
df <- df %>% group_by (Site, Domain, Face, Plate) %>% summarise(freq=n())  #calculate richness by group
df$Method <- "Photo-analysis" #add variable
df$Trait <- "Sessile" #add variable
df$Site<- factor(df$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))
df$Face[df$Face == "B"] <- "Bottom_face" #rename vars
df$Face[df$Face == "T"] <- "Top_face"    #rename vars

# plot site
df$var <- "Site"
p5 <- ggplot (df, aes(x=Site, y=freq, fill=Site))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() +  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(axis.text.x = element_text(size=16),
        legend.position="none", 
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.y = element_text (size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "Species richness") 

#plot face plate
df$var <- "Plate face"
p6 <- ggplot (df, aes(x=Face, y=freq, fill=Face))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#999999"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() +  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(axis.text.x = element_text(size=16),
        legend.position="none", 
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.y = element_text (size=16),
        legend.text =element_text(size=16),
        strip.text = element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "Species richness") 

#plot habitat
df$var <- "Habitat"
p7 <- ggplot (df, aes(x=Domain, y=freq, fill=Domain))  + 
  geom_violin(trim=FALSE, width=1, fill="lightgray") + 
  geom_boxplot(width=0.08, outlier.size = 1) + 
  #geom_jitter(size=3, width=0.06)+
  scale_fill_manual(values=c("#00AFBB","#E7B800"))+
  facet_grid(~var+Method, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(size=16),
        legend.position="none",
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.y = element_text (size=16),
        legend.text=element_text(size=16),
        strip.text = element_text(size=16),
        legend.title=element_blank()) + labs(x ="", y = "Species richnes")


#arrange plots
plot_grid(p_a,p_b,p1,p2,p3,p4,p5,p6,p7, labels = c('a', 'b', "c", "d", "e", "f", "g", "h", "i"), label_size = 22, ncol=2)

#save
ggsave(".././FIGURES/F4.png",  width=17, height=16, dpi=300)

 

