#############################################################################################################
##   Plot ARMS data from 4 datasets.
##   Composition barplots for ARMS identified with genetic and photo-analysis techniques
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################

#load libraries
library (ggplot2)
library (patchwork) 
library (pals)

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")
df_gen$Site  <- factor(df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))
df_photo$Site<- factor(df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1. Plot composition -----
# 1.1 genetic pelagic and benthic phylum
p1 <- ggplot(data=df_gen, aes(x=Size, y=Abundance, fill=Phylum)) + 
  geom_bar(stat='identity',position='fill')+ theme_bw() + 
  scale_fill_manual(values=as.vector(glasbey()))+
  facet_grid(Method~Domain+Site, scales = "free_x")  + 
  theme(legend.text = element_text(size=14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.title=element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=14),
        axis.title.x = element_text(size=14))

# 1.2 genetic pelagic and benthic class
p2 <- ggplot(data=df_gen, aes(x=Size, y=Abundance, fill=Class)) + 
  geom_bar(stat='identity',position='fill')+ 
  facet_grid(Method~Domain+Site, scales = "free_x") +
  scale_fill_manual(values=as.vector(glasbey(32))) + 
  labs(x ="", y = "Relative abundance")+theme_bw() + 
  theme(legend.text = element_text(size=14),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.title=element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=14),
        axis.title.x = element_text(size=14))

#1.3 photo pelagic and benthic phylum
p3 <- ggplot(data=df_photo, aes(x=Trait, y=Cover, fill=Phylum)) + 
  geom_bar(stat='identity',position='fill')+ 
  facet_grid(Method~Domain+Site, scales = "free_x") +
  scale_fill_manual(values=as.vector(glasbey(32))) + 
  labs(x ="", y = "Relative cover") +theme_bw() +
  theme(legend.text = element_text(size=14),
        axis.text = element_text(size = 14),
        legend.title=element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=14),
        axis.title.x = element_text(size=14))

#1.4 photo pelagic and benthic class
p4 <- ggplot(data=df_photo, aes(x=Trait, y=Cover, fill=Class)) + 
  geom_bar(stat='identity',position='fill')+ 
  facet_grid(Method~Domain+Site, scales = "free_x") +
  scale_fill_manual(values=as.vector(glasbey(32))) + 
  labs(x ="", y = "Relative cover")+theme_bw() +
  theme(legend.text = element_text(size=14),
        axis.text = element_text(size = 14),
        legend.title=element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=14),
        axis.title.x = element_text(size=14))

p1/p2/p3/p4+ plot_annotation(tag_levels = 'a') + plot_layout (ncol=1, nrow=4) & theme(plot.tag = element_text(size = 28))

ggsave(".././FIGURES/F3.png",  width=12, height=16, dpi=300)


