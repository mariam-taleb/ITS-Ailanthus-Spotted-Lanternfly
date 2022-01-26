
# Establish environments and libraries -----------------------------------------------
##Libraries
library(vegan)
library(ggplot2)
library(viridis)
library(hillR)
library(readr)
library(dplyr)
library(indicspecies)
library (RColorBrewer)
library (wesanderson)
library(BiodiversityR)
##Set working directory
setwd('C:/Users/Mariam/OneDrive/Documents/Diss. 2017 Data/SM 17/SM17_Pubs_AIL17/Data processing/90%')

##Import basic data
otu_matrix <- read_csv("otu.matrix.csv", 
              col_types = cols(smobs = col_logical()))
otu_matrixfull <- read_csv("otu.matrixfull.csv", 
                       col_types = cols(smobs = col_logical()))
guild_matrix <- read_csv("guild.matrix.csv", 
              col_types = cols(smobs = col_logical()))
tax_frame <- read.csv("tax.frame.csv"), 
              col_types = cols(week = col_number(),smobs = col_logical()))
guild_frame <- read_csv("guild.frame.csv", 
              col_types = cols(week = col_number(),smobs = col_logical()))
gen_matrix <- read_csv("gen.matrix.csv", 
              col_types = cols(week = col_number(),smobs = col_logical()))
gen_order <- c('Trichomerium', 'Strelitziana', 'Sarcinomyces', 'Rhinocladiella', 'Bradymyces', 'Diplodia', 'Dothiorella', 'Alternaria', 'Camarographium', 'Chaetosphaeronema', 'Coniothyrium', 'Curvularia', 'Didymella', 'Epicoccum',   'Leptosphaerulina', 'Lophiostoma', 'Neoascochyta', 'Neocucurbitaria', 'Neosetophoma', 'Neovaginatispora', 'Paraconiothyrium', 'Parapyrenochaeta',   'Phaeosphaeria', 'Phaeosphaeriopsis', 'Phoma', 'Pithomyces', 'Pleospora', 'Sclerostagonospora', 'Setophaeosphaeria', 'Sigarispora', 'Stagonospora',  'Stagonosporopsis', 'Arthrocatena', 'Cladosporium', 'Rachicladosporium', 'Aureobasidium', 'Dothiora', 'Hormonema', 'Leptospora', 'Pringsheimia', 'Sydowia', 'Capnobotryella', 'Cercospora', 'Devriesia', 'Mycosphaerella', 'Neodevriesia', 'Pseudoveronaea', 'Sphaerulina', 'Apiosporina', 'Phyllactinia', 'Articulospora', 'Calloria', 'Phialocephala', 'Phomopsis', 'Plectosphaerella', 'Seiridium', 'Phaeococcomyces', 'Tricellula', 'Taphrina', 'Crepidotus', 'Peniophora', 'Yuchengia', 'Filobasidium', 'Bullera', 'Bulleribasidium', 'Cryptococcus', 'Derxomyces', 'Dimennazyma', 'Dioszegia', 'Fonsecazyma', 'Genolevuria', 'Hannaella', 'Papiliotrema', 'Vishniacozyma', 'Kondoa', 'Buckleyzyma', 'Bannoa', 'Erythrobasidium', 'Symmetrospora', 'Puccinia', 'unidentified')
guild_order<- c('Endophyte- Sooty Mold', 'Endophyte- Undefined Saprobe', 'Litter Saprobe', 'Undefined Saprobe', 'Wood Saprobe', 'Wood Saprobe- Litter Saprobe- Plant Pathogen', 'Litter Saprobe- Animal Pathogen', 'Undefined Saprobe- Plant Pathogen', 'Wood Saprobe- Plant Pathogen', 'Undefined Saprobe- Plant Pathogen- Animal Pathogen',  'Animal Parasite', 'Mycoparasite', 'Plant Pathogen', 'Endophyte- Plant Pathogen- Lichen Parasite', 'Endophyte', 'Lichen', 'Endophyte- Undefined Saprobe- Plant Pathogen- Lichen Parasite', 'Endophyte- Wood Saprobe- Plant Pathogen', 'Endophyte- Wood Saprobe- Plant Pathogen- Fungal Parasite- Lichen Parasite') 
##make data frame Diversity indices
df.div<-data.frame(otu_matrix[,c(1:8)]) #add sample info to div matrix
df.div$shannon <- diversity(otu_matrix[,-c(1:8)], index = "shannon") #shannon
df.div$rich <- hill_taxa(otu_matrix[,-c(1:8)], q=0) #Hill richness
df.div$even <- hill_taxa(otu_matrix[,-c(1:8)], q=2) #Hill evenness

##subset data: include only weeks 4-9
df.div.49<- filter(df.div, week %in% c('4', '5', '6', '7', '8', '9'))
tax_frame49<- filter(tax_frame, week %in% c('4', '5', '6', '7', '8', '9'))
guild_frame49<- filter(guild_frame, week %in% c('4', '5', '6', '7', '8', '9'))
otu_matrix49<- filter(otu_matrix, week %in% c('4', '5', '6', '7', '8', '9'))
otu_matrixfull49<- filter(otu_matrixfull, week %in% c('4', '5', '6', '7', '8', '9'))
gen_matrix49<- filter(gen_matrix, week %in% c('4', '5', '6', '7', '8', '9'))
guild_matrix49<- filter(guild_matrix, week %in% c('4', '5', '6', '7', '8', '9'))

#subset data: all except g=unidentified
tax_frame49conf<- subset(tax_frame49, tax_frame49$g != "unidentified")

##subset only top most abundant genera (by seqct) in weeks 4-9
tax_frame49top<- filter(tax_frame49, g %in% 
                        c('Aureobasidium', 'Cryptococcus',	'Dothiorella',	
                          'Filobasidium',	'Genolevuria',	'Leptosphaerulina',	
                          'Neosetophoma',	'Phoma',	'Phyllactinia',	
                            'Trichomerium'))
guild_frame49top<- filter(guild_frame49, g %in% 
                        c('Aureobasidium', 'Cryptococcus',	'Dothiorella',	
                          'Filobasidium',	'Genolevuria',	'Leptosphaerulina',	
                          'Neosetophoma',	'Phoma',	'Phyllactinia',
                          'Trichomerium'))


# Alpha diversity ---------------------------------------------------------

##normality and transformations
shapiro.test(df.div.49$rich)
histogram(df.div.49$rich)
scatter.smooth(df.div.49$week, df.div.49$rich)


##tests of significance [replace factors as needed]
aov.rich<- aov(rich ~ week, data= df.div.49)
summary(aov.rich)
richtukey<-TukeyHSD(aov.rich)
richtukey
plot(richtukey, las=1, col="#00846b")
kruskal.test(rich ~ week, data= df.div.49)
t.test(df.div.49$rich ~ df.div.49$smobs)

##boxplot alpha div indices [replace factors as needed]
library(ggplot2)
ggplot(data=df.div.49, aes(x=smobs, y=rich, fill=smobs))+
    stat_boxplot(geom = "errorbar", width=0.25, colour = "black") +
  geom_boxplot(width=0.1)+
  geom_jitter(width=0)+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="none")+
  ggtitle("rich-smobs")



##Regressions alpha div by time and weather [replace factors as needed]
#model and significance
linearmodel<-lm(rh ~ even, data=df.div.49)
print(linearmodel)
summary(linearmodel)

#linear regression plots alpha diversity
cor(df.div.49$rich, df.div.49$week)
ggplot(data=df.div.49, aes(x=week, y=rich))+
  geom_point(aes(x=week, y=rich))+
  geom_smooth(method=lm)+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="none")+
  ggtitle("rich-week")


# Relative Abundances-stacked bars *return -----------------------------------------------------

## Genus abundance proportions by week (all genera)
{ggplot(data=tax_frame49, aes(x=factor(week), y=seqct, fill=factor(g, levels = gen_order)))+ 
  geom_bar(position="fill", stat= "identity")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_viridis(discrete=TRUE)+   
  ggtitle("Abundance Proportions by g")+
  xlab("Week")+
  ylab("Abundance Proportion")}


relabunsmobs<- ggplot(data=tax_frame49, aes(x=factor(smobs), y=seqct, fill=factor(g, levels = gen_order)))+ 
  geom_bar(position="fill", stat= "identity")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_viridis(discrete=TRUE)+   
  ggtitle("Abundance Proportions by g")+
  xlab("SM detections")+
  ylab("Abundance Proportion")

ggsave(relabunsmobs, plot= last_plot(), device = "pdf")

## Genus abundance proportions by week (excluding unidentified)
{ggplot(data=tax_frame49conf, aes(x=factor(week), y=seqct, fill=factor(g)))+ 
  geom_bar(position="fill", stat= "identity")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_viridis(discrete=TRUE)+   
  ggtitle("Abundance Proportions by g--unidentified removed")+
  xlab("Week")+
  ylab("Abundance Proportion")}

##Abundance proportions top 10 genera by week
{ggplot(data=guild_frame49, aes(x=factor(week), y=seqct, fill=factor(f)))+ 
  geom_bar(position="fill", stat= "identity")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_viridis(discrete=TRUE)+   
  ggtitle("Abundance Proportions by p")+
  xlab("Week")+
  ylab("Abundance Proportion")}

########### all tax, stacked by week 4-9: MEAN ##not completed--use stat_sum


##Functional group abundance proportions by week
{ggplot(data=guild_frame49, aes(x=factor(week), y=seqct, fill=factor(func, levels = guild_order)))+ 
  geom_bar(position="fill", stat= "identity")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  scale_color_viridis(discrete=TRUE)+
  ggtitle("Abundance Proportions by Functional group by Week")+
  xlab("Week")+
  ylab("Abundance Proportion")}


#Relative abundances linear: week, weather -----------------------------------------------------------------

##Line graph top gen
{ggplot(data=tax_frame49top, aes(week,seqct))+
  geom_line(aes(color = g), size= 1, stat= "summary", fun = mean)+
  geom_point(aes(color = g), size= 2, stat= "summary", fun = mean)+
  theme_classic()+
  #scale_color_brewer(palette = "Set3")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  ggtitle("Mean abundance top 10 genera by week")}

#Line Graph Guilds
{ggplot(data=guild_frame49, aes(week,seqct))+
  geom_line(aes(color = func), size= 1, stat= "summary", fun = mean)+
  geom_point(aes(color = func ), size= 2, stat= "summary", fun = mean)+
  theme_classic()+
  #scale_color_brewer(palette = "Set3")+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  ggtitle("Mean abundance Functional Group by week")}

# Beta diversity by week, weather: absolute abundance ----------------------------------------------------------

##Calculate distance matrix (Bray)
otu49_distmat<- vegdist (otu_matrix49 [,-c(1:8)], method = "bray")
otu49_DistanceMatrix<- as.matrix (gen49_distmat, labels= T)

##Run NMDS
otu49_NMS<- metaMDS(otu_matrix49 [,-c(1:8)], distance = "bray")
stress_otu49<- otu49_NMS$stress

#goodness of fit and shepards
goodness(otu49_NMS)
stressplot(otu49_NMS)

##Create xy data frame for ggplot2 use
{otu49_xy<- data.frame (otu49_NMS$points)
otu49_xy$sample<- otu_matrix49$sample
otu49_xy$week<- otu_matrix49$week
otu49_xy$plant<- otu_matrix49$plant
otu49_xy$rh<- otu_matrix49$rh
otu49_xy$prcp<- otu_matrix49$prcp
otu49_xy$tempmax<- otu_matrix49$tempmax
otu49_xy$tempavg<- otu_matrix49$tempavg
otu49_xy$smobs<- otu_matrix49$smobs}

##check of dispersion
otu49_disper<- betadisper (otu49_distmat, otu_matrix49$week, bias.adjust = TRUE)
permutest(otu49_disper, pairwise= TRUE, permutations = 99)

#dispersion data frame
disper_dist<-data.frame(otu_matrix49 [c(1:8)])
disper_dist$distance<-otu49_disper$distances   

#Smobs~dispersions
t.test(disper_dist$distance, disper_dist$smobs)
ggplot(data=disper_dist, aes(x=smobs, y=distance, fill=smobs))+
  stat_boxplot(geom = "errorbar", width=0.25, colour = "black") +
  geom_boxplot(width=0.1)+
  geom_jitter(width=0)+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="none")+
  ggtitle("Dispersals by SM Detections")

#week~ dispersions
print(lm(distance ~ week, data=disper_dist))
summary(lm(distance ~ week, data=disper_dist))
ggplot(data=disper_dist, aes(x=week, y=distance))+
  geom_point(aes(x=week, y=distance))+
  geom_smooth(method=lm)+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="none")+
  ggtitle("Dispersals over Time")



##PERMANOVA (adonis2)
otu49_adonis<- adonis2(formula= otu_matrix49[,-c(1:8)]  ~ otu_matrix49$week, 
              method="bray", perm= 999)
print(otu49_adonis)

##envfit
otu49_env<- data.frame(otu_matrix49 [c(1:8)])
otu49_envfit<- envfit(otu49_NMS, otu49_env[,-c(3)] , permutations = 999)
print(otu49_envfit)
otu49_vector_coord <- as.data.frame(scores(otu49_envfit, "vectors"))* 
  ordiArrowMul(otu49_envfit)
otu49_factor_coord <- as.data.frame(scores(otu49_envfit, "factors"))* 
  ordiArrowMul(otu49_envfit)

##Plot ordination:ggplot  
{ggplot(otu49_xy, aes(x=MDS1, y=MDS2))+
  geom_point(size=4, shape= "triangle", alpha= 0.6, 
    aes(group=factor(week), color=factor(week)))+ #plot all samples
  geom_text(data = otu49_xy, aes(x = MDS1, y = MDS2+0.02), 
    label = otu49_xy$sample, size= 2, alpha= 1, 
    colour = "black") + #label all samples
  scale_color_viridis(discrete=TRUE)+ #Viridis palette for samples
  stat_ellipse(aes(x=MDS1, y=MDS2, group=factor(smobs), linetype=factor(smobs)), 
     size=1, data = NULL, geom = "path", position = "identity", type = "t", 
     level = 0.95, segments = 51)+ #ellipses by factor
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
      data = otu49_vector_coord[1,], size=1, alpha= 1, colour = "black", 
      arrow = arrow(length=unit(6, "points"))) + #envfit vector(p<0.05)
  geom_text(data = otu49_vector_coord[1,], aes(x=NMDS1, y=NMDS2-0.01),size=4,
      colour="black",label=row.names(otu49_vector_coord[1,])) + #label vector
  geom_point(data = otu49_factor_coord [6:7,], aes(x = NMDS1, y = NMDS2), 
      shape = "circle", size = 6, alpha = .3, colour = "black") + #centroids SM
  geom_text(data =otu49_factor_coord[6:7,],aes(x=NMDS1, y=NMDS2-0.02), label= 
      row.names(otu49_factor_coord[6:7,]),size=3,colour="black")+ #label cent SM
  #geom_point(data =otu49_factor_coord [1:5,], aes(x = NMDS1, y = NMDS2), 
      #shape = "square", size = 3, alpha = .5, colour = "black") + #cent plant
  #geom_text(data = otu49_factor_coord[1:5,], aes(x = NMDS1, y = NMDS2+0.01), 
      #label = row.names(otu49_factor_coord[1:5,]), size= 2, alpha= 0.7, 
      #colour = "black") + #label centroids plant
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))}
  
# Beta diversity :all OTUs ------------------------------------------------
##Calculate distance matrix (Bray)
otufull49_distmat<- vegdist (otu_matrixfull49 [,-c(1:8)], method = "bray")
otufull49_DistanceMatrix<- as.matrix (gen49full_distmat, labels= T)

##Run NMDS
otufull49_NMS<- metaMDS(otu_matrixfull49 [,-c(1:8)], distance = "bray")
stress_otufull49<- otufull49_NMS$stress

#goodness of fit and shepards
goodness(otufull49_NMS)
stressplot(otufull49_NMS)

##Create xy data frame for ggplot2 use
{otufull49_xy<- data.frame (otufull49_NMS$points)
  otufull49_xy$sample<- otu_matrixfull49$sample
  otufull49_xy$week<- otu_matrixfull49$week
  otufull49_xy$plant<- otu_matrixfull49$plant
  otufull49_xy$rh<- otu_matrixfull49$rh
  otufull49_xy$prcp<- otu_matrixfull49$prcp
  otufull49_xy$tempmax<- otu_matrixfull49$tempmax
  otufull49_xy$tempavg<- otu_matrixfull49$tempavg
  otufull49_xy$smobs<- otu_matrixfull49$smobs 
  }

##check of dispersion
otufull49_disper<- betadisper (otufull49_distmat, otu_matrixfull49$week)
permutest(otufull49_disper, pairwise= TRUE, permutations = 99)
otufull49_disper

##PERMANOVA (adonis2)
otufull49_adonis<- adonis2(formula= otu_matrixfull49[,-c(1:8)]  ~ otu_matrixfull49$prcp, 
                       method="bray", perm= 999)
print(otufull49_adonis)

##envfit
otufull49_env<- data.frame(otu_matrixfull49 [c(1:8)])
otufull49_envfit<- envfit(otufull49_NMS, otufull49_env[,-c(3)] , permutations = 999)
print(otufull49_envfit)
otufull49_vector_coord <- as.data.frame(scores(otufull49_envfit, "vectors"))* 
  ordiArrowMul(otufull49_envfit)
otufull49_factor_coord <- as.data.frame(scores(otufull49_envfit, "factors"))* 
  ordiArrowMul(otufull49_envfit)

##Plot ordination:ggplot  
{ggplot(otufull49_xy, aes(x=MDS1, y=MDS2))+
    geom_point(size=4, shape= "triangle", alpha= 0.6, 
               aes(group=factor(week), color=factor(week)))+ #plot all samples
    geom_text(data = otufull49_xy, aes(x = MDS1, y = MDS2+0.02), 
              label = otufull49_xy$sample, size= 2, alpha= 1, 
              colour = "black") + #label all samples
    scale_color_viridis(discrete=TRUE)+ #Viridis palette for samples
    stat_ellipse(aes(x=MDS1, y=MDS2, group=factor(smobs), linetype=factor(smobs)), 
                 size=1, data = NULL, geom = "path", position = "identity", type = "t", 
                 level = 0.95, segments = 51)+ #ellipses by factor
    #geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 #data = otufull49_vector_coord[1,], size=1, alpha= 1, colour = "black", 
                 #arrow = arrow(length=unit(6, "points"))) + #envfit vector(p<0.05)
    #geom_text(data = otufull49_vector_coord[1,], aes(x=NMDS1, y=NMDS2-0.01),size=4,
             # colour="black",label=row.names(otu49_vector_coord[1,])) + #label vector
   # geom_point(data = otufull49_factor_coord [6:7,], aes(x = NMDS1, y = NMDS2), 
    #           shape = "circle", size = 6, alpha = .3, colour = "black") + #centroids SM
  #  geom_text(data =otufull49_factor_coord[6:7,],aes(x=NMDS1, y=NMDS2-0.02), label= 
    #            row.names(otu49_factor_coord[6:7,]),size=3,colour="black")+ #label cent SM
    #geom_point(data =otu49_factor_coord [1:5,], aes(x = NMDS1, y = NMDS2), 
    #shape = "square", size = 3, alpha = .5, colour = "black") + #cent plant
    #geom_text(data = otu49_factor_coord[1:5,], aes(x = NMDS1, y = NMDS2+0.01), 
    #label = row.names(otu49_factor_coord[1:5,]), size= 2, alpha= 0.7, 
    #colour = "black") + #label centroids plant
    theme_classic()+
    theme(panel.background = element_rect(fill = "transparent", colour = NA),  
          plot.background = element_rect(fill = "transparent", colour = NA))
  }

# Beta diversity by week, weather: relative abundance (**VOID:See notes) ---------------------
      ###**BC is decently tolerant of the kinds of distribution in OTU tables, 
      ###*but relative abundance can make for better analyses in some cases.
      ###*In this instance, stress was too low, and indicated relative abundance 
      ###*is a poor choice for this data set 

# community dissimilary tests ---------------------------------------------

#model and significance
BCline<-lm( ~ even, data=df.div.49)
print(linearmodel)
summary(linearmodel)

#linear regression plots alpha diversity
cor(df.div.49$rich, df.div.49$week)
ggplot(data=df.div.49, aes(x=week, y=rich))+
  geom_point(aes(x=week, y=rich))+
  geom_smooth(method=lm)+
  theme_classic()+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="none")+
  ggtitle("rich-week")
# Indicator Species -------------------------------------------------------
### Indval to examine relation of genera to site group: smobs
gen_smobs<-multipatt(gen_matrix49 [,-(1:8)], gen_matrix49$smobs, 
    func = "r.g", control = how(nperm=9999))
summary (gen_smobs)

### Indval to examine relation of functional group to site group: smobs
guild_smobs<-multipatt(guild_matrix49 [,-(1:8)], guild_matrix49$smobs, 
    func = "r.g", control = how(nperm=9999))
summary (guild_smobs)

### Indval to examine relation of genera to site group: week
gen_week<-multipatt(gen_matrix49 [,-(1:8)], gen_matrix49$week, 
                     func = "r.g", control = how(nperm=9999))
summary (gen_week)

### Indval to examine relation of functional group to site group: week
guild_week<-multipatt(guild_matrix49 [,-(1:8)], guild_matrix49$week, 
                       func = "r.g", control = how(nperm=9999))
summary (guild_week)

# Co-occurence models -----------------------------------------------------
### examine competition between genera
###examine competition between funcitonal groups
