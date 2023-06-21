#############################################################################################################
##   ANOSIM in pelagic and benthic habitats
##   Project: IM20-UrbanKlima.
##   Author: E. Villarino, evillarino@azti.es
##############################################################################################################


## 0. load libraries
library (vegan) 
library (tidyverse)
library (viridis)
library (cowplot)
library (ggplot2)

#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# 0. Load dataframes previously saved in script 1 
load("./dat_arms.RData")
df_gen$Site  <- factor(df_gen$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))
df_photo$Site<- factor(df_photo$Site, levels = c("Mendexa-pelagic","Lekeitio","Zumaia", "Pasaia"))

# 1. manipulate data photo-analysis 
df_photo_site_face <- df_photo  %>% 
  group_by(Site,Species,Face,Plate,Domain) %>% 
  summarize(Cover = mean(Cover, na.rm = TRUE)) %>% 
  spread(Species,Cover) %>% 
  select (-Unidentified)

# 1.1 Calculate Bray-Curtis photo-analysis
df_bc_photo_site_face <- vegdist (df_photo_site_face[5:33], "bray", na.rm = T)  # Bray-Curtis.

# 2. manipulate data genetic
df_gen_size_site <-  df_gen   %>% 
  group_by(Site,SWARM,Size,Domain, Trait) %>% 
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  spread(SWARM,Abundance)

# 2.1 Calculate Bray-Curtis genetic
df_bc_gen_size_site <- vegdist (df_gen_size_site[5:560], "bray", na.rm = T)  # Bray-Curtis.

# 3. Calculate ANOSIM photo-analysis -----
# ("The ANOSIM statistic "R" compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups)
test1 <- anosim(df_bc_photo_site_face, df_photo_site_face$Site,  permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.59 # significance: 0.001
test2 <- anosim(df_bc_photo_site_face, df_photo_site_face$Domain,permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.87 # significance: 0.001
test3 <- anosim(df_bc_photo_site_face, df_photo_site_face$Plate, permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=-0.06 # significance: 0.895
test4 <- anosim(df_bc_photo_site_face, df_photo_site_face$Face,  permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.01 # significance: 0.286

# 3.1 Dissimilarity rank photo
df1 <- cbind.data.frame(Dis_rank=test1$dis.rank,Site=test1$class.vec)
df2 <- cbind.data.frame(Dis_rank=test2$dis.rank,Domain=test2$class.vec)
df3 <- cbind.data.frame(Dis_rank=test3$dis.rank,Plate=test3$class.vec)
df4 <- cbind.data.frame(Dis_rank=test4$dis.rank,Face=test4$class.vec)

# 3.2 Add variable
df1$Method <- "Photo-analysis"
df2$Method <- "Photo-analysis"
df3$Method <- "Photo-analysis"
df4$Method <- "Photo-analysis"

# 3.3 Rename variable
levels(df4$Face)[levels(df4$Face)=='T'] <- 'Top'
levels(df4$Face)[levels(df4$Face)=='B'] <- 'Bottom'

# 4. Calculate ANOSIM genetic-----
test5 <- anosim(df_bc_gen_size_site, df_gen_size_site$Site,   permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.66 # significance: 0.001
test6 <- anosim(df_bc_gen_size_site, df_gen_size_site$Domain, permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.82 # significance: 0.002
test7 <- anosim(df_bc_gen_size_site, df_gen_size_site$Size,   permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.12 # significance: 0.162
test8 <- anosim(df_bc_gen_size_site, df_gen_size_site$Trait,  permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))    # R=0.01 # significance: 0.43

# 4.1 Dissimilarity rank genetic
df5 <- cbind.data.frame(Dis_rank=test5$dis.rank,Site=test5$class.vec)
df6 <- cbind.data.frame(Dis_rank=test6$dis.rank,Domain=test6$class.vec)
df7 <- cbind.data.frame(Dis_rank=test7$dis.rank,Size=test7$class.vec)
df8 <- cbind.data.frame(Dis_rank=test8$dis.rank,Trait=test8$class.vec)

# 4.2. Add variable
df5$Method <- "Metabarcoding"
df6$Method <- "Metabarcoding"
df7$Method <- "Metabarcoding"
df8$Method <- "Metabarcoding"

#5. plot photo-analysis-----

# add labels to plot
dat_stats1 <- data.frame (Label=c("R=0.59  p=0.001"),Site =rep(c("Site"),times=1),x=2.5,y=350)
dat_stats2 <- data.frame (Label=c("R=0.87  p=0.001"),Domain =rep(c("Domain"),times=1),x=2.5,y=350)
dat_stats3 <- data.frame (Label=c("R=-0.06 p>0.05"),Plate =rep(c("Plate"),times=1),x=2.5,y=350)
dat_stats4 <- data.frame (Label=c("R=0.01  p>0.05"),Face =rep(c("Face"),times=1),x=2.5,y=350)

#site
df1$var <- "Site"
p1 <- ggplot(df1,aes(x=Site, y=Dis_rank, fill=Site)) +
  geom_violin(trim=FALSE, fill="gray", width=1.4) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Site))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats1,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values=c("Between"="dodgerblue","Mendexa-pelagic"="#00AFBB", "Lekeitio"="#E7B800", "Zumaia"="#FC4E07", "Pasaia"="#999999"))+
  scale_color_manual(values=c("Between"="dodgerblue","Mendexa-pelagic"="#00AFBB", "Lekeitio"="#E7B800", "Zumaia"="#FC4E07", "Pasaia"="#999999"))#jitter color

#habitat
df2$var <- "Habitat"
p2 <- ggplot(df2,aes(x=Domain, y=Dis_rank, fill=Domain)) +
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Domain))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats2,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+
  scale_fill_manual(values=c("Between"="dodgerblue","Pelagic"="#00AFBB","Benthic"="#E7B800"))+
  scale_color_manual(values=c("Between"="dodgerblue","Pelagic"="#00AFBB","Benthic"="#E7B800"))#jitter color

#plate number
df3$var <- "Plate number"
p3 <- ggplot(df3,aes(x=Plate, y=Dis_rank, fill=Plate)) +
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Plate))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats3,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+
  scale_fill_manual(values=c("Between"="dodgerblue","1"="#00AFBB", "4"="#E7B800", "8"="#FC4E07"))+
  scale_color_manual(values=c("Between"="dodgerblue","1"="#00AFBB", "4"="#E7B800", "8"="#FC4E07"))#jitter color

#plate face
df4$var <- "Plate face"
p4 <- ggplot(df4,aes(x=Face, y=Dis_rank, fill=Face)) +
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") +  
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Face))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats4,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+
  scale_fill_manual(values=c("Between"="dodgerblue","Bottom"="#00AFBB","Top"="#E7B800"))+
  scale_color_manual(values=c("Between"="dodgerblue","Bottom"="#00AFBB","Top"="#E7B800"))#jitter color

#6. plot genetic ----

#add labels to plot
dat_stats5 <- data.frame (Label=c("R=0.66  p=0.001"),Site =rep(c("Site"),times=1),x=2.5,y=240)
dat_stats6 <- data.frame (Label=c("R=0.82  p=0.002"),Domain =rep(c("Domain"),times=1),x=2.5,y=240)
dat_stats7 <- data.frame (Label=c("R=0.12  p>0.05"),Size =rep(c("Size"),times=1),x=2.5,y=240)
dat_stats8 <- data.frame (Label=c("R=0.01  p>0.05"),Trait =rep(c("Trait"),times=1),x=2.5,y=240)

#site
df5$var <- "Site"
p5 <- ggplot(df5,aes(x=Site, y=Dis_rank, fill=Site)) +
  geom_violin(trim=FALSE, fill="gray", width=1.2) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Site))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats5,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_fill_manual(values=c("Between"="dodgerblue","Mendexa-pelagic"="#00AFBB", "Lekeitio"="#E7B800", "Zumaia"="#FC4E07", "Pasaia"="#999999"))+
  scale_color_manual(values=c("Between"="dodgerblue","Mendexa-pelagic"="#00AFBB", "Lekeitio"="#E7B800", "Zumaia"="#FC4E07", "Pasaia"="#999999"))#jitter color

#habitat
df6$var <- "Habitat"
p6 <- ggplot(df6,aes(x=Domain, y=Dis_rank, fill=Domain)) +
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Domain))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats6,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+
  scale_fill_manual(values=c("Between"="dodgerblue","Pelagic"="#00AFBB","Benthic"="#E7B800"))+
  scale_color_manual(values=c("Between"="dodgerblue","Pelagic"="#00AFBB","Benthic"="#E7B800"))#jitter color

#size
df7$var <- "Size"
p7 <- ggplot(df7,aes(x=Size, y=Dis_rank, fill=Size)) +
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Size))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats7,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+
  scale_fill_manual(values= c("Between"="dodgerblue",">2000"="#00AFBB","500-1000"="#E7B800", "90-500"="#FC4E07", "Sessile"="#999999"))+
  scale_color_manual(values=c("Between"="dodgerblue",">2000"="#00AFBB","500-1000"="#E7B800", "90-500"="#FC4E07", "Sessile"="#999999"))

#trait
df8$var <- "Trait"
p8 <- ggplot(df8,aes(x=Trait, y=Dis_rank, fill=Trait)) +
  geom_violin(trim=FALSE, fill="gray", width=1) +
  geom_boxplot(width=0.07, size=0.6, color="black") + 
  theme_bw() + facet_grid(~var+Method)+
  geom_jitter(size=5, width=0.25, alpha = 0.25, aes(size = 2, color=Trait))+
  xlab("") + ylab("Dissimilarity rank")+ geom_text(data = dat_stats8,aes(x=x, y=y, label = Label),color="black",size=5)+
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size=16), # remove x-axis labels
        axis.title.y = element_text(size=16), # remove y-axis labels
        axis.text.x = element_text (size=14),
        axis.text.y = element_text (size=16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.position = "none")+
  scale_fill_manual(values=c("Between"="dodgerblue","Mobile"="#00AFBB","Sessile"="#E7B800"))+
  scale_color_manual(values=c("Between"="dodgerblue","Mobile"="#00AFBB","Sessile"="#E7B800"))#jitter color

plot_grid(p5,p6,p7,p8, p1,p2,p3,p4, labels = c('a', 'b', "c", "d", "e", "f", "g", "h"), label_size = 22, ncol=2)


ggsave(".././FIGURES/F7.png",  width=10, height=12, dpi=300)




