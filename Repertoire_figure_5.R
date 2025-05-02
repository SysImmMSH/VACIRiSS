##### updated vaciriss repertroire script for paper 1 #######

##LOAD PACKAGES##
library(immunarch)
library(tidyverse)
library(ggplot2)
library(circlize)
library(ggpubr)
#install.packages("rstatix")
library(rstatix)
## working directory 

getwd()

immdata<- repLoad("Sample_report_2/") 


######### B cells #########
bcell<- repFilter(immdata, "by.clonotype", list(V.name= include("IG")), .match="substring")
#repSave(bcell, "Repertoire/")

meta <- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata_repertoire.csv")

meta
meta$Treatment[meta$Treatment == 'B'] <- 'PCV13'
meta$Treatment[meta$Treatment == 'A'] <- 'Placebo'
meta

meta$Status<- paste(meta$Treatment, meta$Timepoint)
meta


colnames(meta)

colnames(meta)
meta

meta$Cohort_time<- paste(meta$Treatment, meta$Timepoint, sep = "_")

bcell$meta
bcell


bcell$meta<- meta
bcell$meta
bcell

### Clonality volume#####
###### When .method is set to "volume" the repExplore calculates the number of unique clonotypes in the input data.#####
clono_vol_b<- repExplore(bcell$data, .method = "volume")
clono_vol_b<- as.data.frame(clono_vol_b)

meta_clone_volume<-merge(x = clono_vol_b, y = meta, by= "Sample", all = TRUE)
meta_clone_volume

write.csv(meta_clone_volume, "BCell_clono_number_volume_T0_T1_T2_anon.csv")

clono_vol_b_df<- merge(clono_vol_b, meta)

write.csv(clono_vol_b_df, "BCell_clono_number_volume_T0_T1_T2_anon.csv")

clono_vol_b_df<-read.csv("BCell_clono_number_volume_T0_T1_T2_anon.csv")
clono_vol_b_df<- clono_vol_b_df[!clono_vol_b_df$Timepoint == "T3", ]
clono_vol_b_df<- clono_vol_b_df[!clono_vol_b_df$Timepoint == "T1", ]


#### outer merge to check missing samples

meta_clone_volume<-merge(x = clono_vol_b_df, y = meta, by= "Sample", all = TRUE)
meta_clone_volume



clono_vol_b_df$Treatment <- factor(clono_vol_b_df$Treatment, levels=c('Placebo', 'PCV13'))

clono_vol_b_df$Timepoint <- as.factor(clono_vol_b_df$Timepoint)

#clono_vol_b_df$Status <- factor(clono_vol_b_df$Status, levels=c('T0 PCV13', 'T1 PCV13', 'T2 PCV13', 'T3 PCV13', 'T0 Placebo', 'T1 Placebo', 'T2 Placebo', 'T3 Placebo'))
clono_vol_b_df$Status <- factor(clono_vol_b_df$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0', 'Placebo T2'))

clono_vol_b_df$Status <- as.factor(clono_vol_b_df$Status)
clono_vol_b_df
colnames(clono_vol_b_df)

my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"),  c("PCV13 T0", "PCV13 T2"))


clono_plot<-ggplot(clono_vol_b_df, aes(x= Status, y= Volume)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("Unique clones")+
  xlab("")+
  ggtitle("Volume")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_clones, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot
ggsave(plot = clono_plot, "Figures/Repertoire/B cell/Paper/Clone_volume_unique.png", dpi = 300,  height = 6,  width = 6)


######## get numbers of clones for each status 


clono_vol_b_df

summary_stats <- clono_vol_b_df %>%
  group_by(Status) %>%
  summarise(
    Mean = mean(Volume),
    Median = median(Volume),
    Q1 = quantile(Volume, 0.25),
    Q3 = quantile(Volume, 0.75)
  )

summary_stats



#####################
### add responder here 

clono_vol_b_df

#responder<- read.csv("Responder/Responder_patientID.csv")
#responder

clono_vol_b_df

#clono_vol_b_df_resp<- merge(clono_vol_b_df, responder, by= "patient_ID")


clono_vol_b_df_resp<- clono_vol_b_df[!is.na(clono_vol_b_df$Responder_final),]
clono_vol_b_df_resp

table(clono_vol_b_df_resp$Responder_final)

clono_vol_b_df_resp
clono_vol_b_df_resp$Timepoint_resp<- paste(clono_vol_b_df_resp$Timepoint, clono_vol_b_df_resp$Responder_final, sep = "_")

my_comps3 <- list(c("T0_Yes", "T0_No"), c("T2_Yes", "T2_No"))


clono_plot<-ggplot(clono_vol_b_df_resp, aes(x= Timepoint_resp, y= Volume)) +
  geom_boxplot(aes(fill= Timepoint_resp), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Timepoint_resp), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("T0_Yes" = "#f58f71",
                              "T2_Yes" = "#f58f71",
                              "T0_No" ="#02dbc7",
                              "T2_No" ="#02dbc7"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("Unique clones")+
  xlab("")+
  ggtitle("Volume")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black")# +
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps3, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot
ggsave(plot = clono_plot, "Figures/Repertoire/B cell/Paper/Clone_volume_unique_responder_Bcell_timepoint.png", dpi = 300,  height = 6,  width = 6)


summary_stats <- clono_vol_b_df_resp %>%
  group_by(Timepoint_resp) %>%
  summarise(
    Mean = mean(Volume),
    Median = median(Volume),
    Q1 = quantile(Volume, 0.25),
    Q3 = quantile(Volume, 0.75)
  )

summary_stats
table(clono_vol_b_df_resp$Timepoint_resp)

####################

####  diversity #######
####### Diversity choa ####

div_chao_b <- repDiversity(bcell$data, "chao1")
div_chao_b
div_chao_b<- as.data.frame(div_chao_b)

div_chao_b

div_chao_b$Sample <- rownames(div_chao_b)
rownames(div_chao_b)<- NULL

div_chao_b<- merge(div_chao_b, meta)

div_chao_b

div_chao_b<- div_chao_b[!div_chao_b$Timepoint == "T1", ]
div_chao_b<- div_chao_b[!div_chao_b$Timepoint == "T3", ]

div_chao_b$Status <- factor(div_chao_b$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0', 'Placebo T2'))


my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"),  c("PCV13 T0", "PCV13 T2"))


clono_plot<-ggplot(div_chao_b, aes(x= Status, y= Estimator)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  #ylab("Clones")+
  xlab("")+
  ggtitle("Chao 1")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_clones, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot

#c("Placebo T0", "Placebo T2"),
ggsave(plot = clono_plot, "Figures/Repertoire/B cell/Paper/Chao_diversity.png", dpi = 300,  height = 6,  width = 6)


summary_stats <- div_chao_b %>%
  group_by(Status) %>%
  summarise(
    Mean = mean(Estimator),
    Median = median(Estimator),
    Q1 = quantile(Estimator, 0.25),
    Q3 = quantile(Estimator, 0.75)
  )

summary_stats


####### Diversity simp ####

div_simp_b <- repDiversity(bcell$data, "inv.simp")
div_simp_b
div_simp_b<- as.data.frame(div_simp_b)

div_simp_b


meta
div_simp_b<- merge(div_simp_b, meta)

div_simp_b

div_simp_b<- div_simp_b[!div_simp_b$Timepoint == "T3", ]
div_simp_b<- div_simp_b[!div_simp_b$Timepoint == "T1", ]



div_simp_b$Status <- factor(div_simp_b$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0', 'Placebo T2'))


my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"))

colnames(div_simp_b)

clono_plot<-ggplot(div_simp_b, aes(x= Status, y= Value)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("1/D")+
  xlab("")+
  ggtitle("Inverse Simpson Index")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_clones, method= "t.test",  label = "p.signif", hide.ns = FALSE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot

#c("Placebo T0", "Placebo T2"),
#ggsave(plot = clono_plot, "Figures/Repertoire/B cell/inver_simp_diversity.png", dpi = 300,  height = 6,  width = 6)
ggsave(plot = clono_plot, "Figures/Repertoire/B cell/Paper/inver_simp_diversity.png", dpi = 300,  height = 6,  width = 6)

summary_stats <- div_simp_b %>%
  group_by(Status) %>%
  summarise(
    Mean = mean(Value),
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

summary_stats

######## responder diversity

div_simp_b_resp<- div_simp_b[!is.na(div_simp_b$Responder_final),]


div_simp_b_resp
div_simp_b_resp$Timepoint_resp<- paste(div_simp_b_resp$Timepoint, div_simp_b_resp$Responder_final, sep = "_")

my_comps3 <- list(c("T0_Yes", "T0_No"), c("T2_Yes", "T2_No"))

clono_plot<-ggplot(div_simp_b_resp, aes(x= Timepoint_resp, y= Value)) +
  geom_boxplot(aes(fill= Timepoint_resp), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Timepoint_resp), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("T0_Yes" = "#f58f71",
                              "T2_Yes" = "#f58f71",
                              "T0_No" ="#02dbc7",
                              "T2_No" ="#02dbc7"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("1/D")+
  xlab("")+
  ggtitle("Inverse Simpson Index")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps3, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot

#c("Placebo T0", "Placebo T2"),
#ggsave(plot = clono_plot, "Figures/Repertoire/B cell/inver_simp_diversity.png", dpi = 300,  height = 6,  width = 6)
ggsave(plot = clono_plot, "Figures/Repertoire/B cell/Paper/inver_simp_diversity_overall_Responder_Bcell_timepoint.png", dpi = 300,  height = 6,  width = 6)

summary_stats <- div_simp_b_resp %>%
  group_by(Timepoint_resp) %>%
  summarise(
    Mean = mean(Value),
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

summary_stats


##########################
###### Rep clonality ####
##########################


# Set "homeo" to analyse relative abundance (also known as clonal space homeostasis),
#which is defined as the proportion of repertoire occupied by clonal groups with specific abundances..###

imm_hom <- repClonality(bcell$data, .method = "homeo")

imm_hom <- as.data.frame(imm_hom)
imm_hom

imm_hom$Sample <- rownames(imm_hom)
rownames(imm_hom)<- NULL

imm_hom
imm_hom_df<- merge(meta, imm_hom)
imm_hom_df


library(tidyr)

imm_hom_df_piv<-imm_hom_df %>% pivot_longer(cols=c('Rare (0 < X <= 1e-05)', 'Small (1e-05 < X <= 1e-04)', 'Medium (1e-04 < X <= 0.001)', 'Large (0.001 < X <= 0.01)', 'Hyperexpanded (0.01 < X <= 1)'),
                                            names_to='Clone_size',
                                            values_to='Values')


imm_hom_df_piv
imm_hom_df_piv<- imm_hom_df_piv[!imm_hom_df_piv$Timepoint == "T1", ]
imm_hom_df_piv<- imm_hom_df_piv[!imm_hom_df_piv$Timepoint == "T3", ]


imm_hom_df_piv$Status <- factor(imm_hom_df_piv$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0',  'Placebo T2'))


my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"), c("PCV13 T0", "PCV13 T2"))

colnames(imm_hom_df_piv)



imm_hom_df_piv$Clone_size <- factor(imm_hom_df_piv$Clone_size, levels=c('Rare (0 < X <= 1e-05)', 'Small (1e-05 < X <= 1e-04)', 'Medium (1e-04 < X <= 0.001)', 'Large (0.001 < X <= 0.01)', 'Hyperexpanded (0.01 < X <= 1)'))


mycomps_1<- list( c("T0 PCV13", "T1 PCV13"), c("T0 PCV13", "T2 PCV13"), c("T0 PCV13", "T3 PCV13"), c("T0 Placebo", "T1 Placebo"), c("T0 Placebo", "T2 Placebo"), c("T0 Placebo", "T3 Placebo"))

imm_hom_df_piv$Status <- as.factor(imm_hom_df_piv$Status)
imm_hom_df_piv$Treatment <- factor(imm_hom_df_piv$Treatment, levels=c('PCV13','Placebo'))


imm_hom_df_piv$Clone_size_status<- paste(imm_hom_df_piv$Clone_size, imm_hom_df_piv$Cohort_time, sep = "_")

my_comps_clones <- list( c("Placebo T0", "Placebo T2"), c("PCV13 T0", "PCV13 T2"))


# Split data by species
iris_split <- split(imm_hom_df$`Medium (1e-04 < X <= 0.001)`, imm_hom_df$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if greater than 0.05 then it is not normal 



####### with facet ####### 
clono_plot<- ggplot(imm_hom_df_piv, aes(x= Clone_size, y= Values, fill=Status))+ 
  geom_boxplot(aes(fill= Status), outlier.shape = NA, shape =21, alpha = 0.4)+
  #geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 1) +
  geom_point(aes(fill=Status), shape= 21, position = position_jitterdodge(), size = 3)+
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("Percentage")+
  xlab("")+
  ggtitle("Clonal distribution")+
  theme(plot.title = element_text(hjust=0.5)) +
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  #geom_vline(xintercept = 3.5, linetype= "dotted", color = "black") +
  facet_wrap(vars(Treatment))+
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps2, method= "t.test",  label = "p.signif", hide.ns = TRUE)
  geom_signif(comparisons = my_comps2, map_signif_level = TRUE, test = "wilcox.test") 

clono_plot
w <-clono_plot + geom_pwc(aes(group = Status), method = "wilcox_test", label = "p.signif", hide.ns = TRUE)
w

ggsave("Clonality_homeo_stats_nonfacet_wilcox.png", w, path ="Figures/Repertoire/B cell/Paper/",  width = 20, height= 15, dpi = 300 )



##########################################
##########################################
############# T cells ####################
##########################################
##########################################

tcell<- repFilter(immdata, "by.clonotype", list(V.name= include("TR")), .match="substring")

meta <- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata_repertoire.csv")

#meta
meta$Treatment[meta$Treatment == 'B'] <- 'PCV13'
meta$Treatment[meta$Treatment == 'A'] <- 'Placebo'
meta


meta$Status<- paste(meta$Treatment, meta$Timepoint)
#meta

tcell$meta<- meta
tcell$meta
tcell

#repSave(bcell)

#immdata$meta<- meta
### Clonality volume#####
###### When .method is set to "volume" the repExplore calculates the number of unique clonotypes in the input data.#####
clono_vol_t<- repExplore(tcell$data, .method = "volume")
clono_vol_t<- as.data.frame(clono_vol_t)

meta_clone_volume_t<-merge(x = clono_vol_t, y = meta, by= "Sample", all = TRUE)
meta_clone_volume_t

write.csv(meta_clone_volume_t, "TCell_clono_number_volume_T0_T1_T2_anon.csv")

clono_vol_t_df<- merge(clono_vol_t, meta)

write.csv(clono_vol_t_df, "TCell_clono_number_volume_T0_T1_T2_anon.csv")

clono_vol_t_df<-read.csv("TCell_clono_number_volume_T0_T1_T2_anon.csv")
clono_vol_t_df<- clono_vol_t_df[!clono_vol_t_df$Timepoint == "T3", ]
clono_vol_t_df<- clono_vol_t_df[!clono_vol_t_df$Timepoint == "T1", ]


#### outer merge to check missing samples

meta_clone_volume<-merge(x = clono_vol_t_df, y = meta, by= "Sample", all = TRUE)
meta_clone_volume

########

clono_vol_t_df$Treatment <- factor(clono_vol_t_df$Treatment, levels=c('PCV13', 'Placebo'))

clono_vol_t_df$Timepoint <- as.factor(clono_vol_t_df$Timepoint)

#clono_vol_b_df$Status <- factor(clono_vol_b_df$Status, levels=c('T0 PCV13', 'T1 PCV13', 'T2 PCV13', 'T3 PCV13', 'T0 Placebo', 'T1 Placebo', 'T2 Placebo', 'T3 Placebo'))
clono_vol_t_df$Status <- factor(clono_vol_t_df$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0', 'Placebo T2'))

clono_vol_t_df$Status <- as.factor(clono_vol_t_df$Status)
clono_vol_t_df
colnames(clono_vol_t_df)

clono_vol_t_df

my_comps_clones <- list( c("Placebo T0", "Placebo T2"),  c("PCV13 T0", "PCV13 T2"))


clono_plot<-ggplot(clono_vol_t_df, aes(x= Status, y= Volume)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("Unique clones")+
  xlab("")+
  ggtitle("Volume")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_clones, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot


ggsave(plot = clono_plot, "Figures/Repertoire/T cell/Paper/Clone_volume_unique.png", dpi = 300,  height = 6,  width = 6)



######## get numbers of clones for each status 


clono_vol_t_df


summary_stats_T <- clono_vol_t_df %>%
  group_by(Status) %>%
  summarise(
    Mean = mean(Volume),
    Median = median(Volume),
    Q1 = quantile(Volume, 0.25),
    Q3 = quantile(Volume, 0.75)
  )

summary_stats_T
### add responder here 
clono_vol_t_df

clono_vol_t_df_resp<- clono_vol_t_df[!is.na(clono_vol_t_df$Responder_final),]


my_comps2 <- list(c("Yes", "No"))


clono_vol_t_df_resp$Timepoint_resp<- paste(clono_vol_t_df_resp$Timepoint, clono_vol_t_df_resp$Responder_final, sep = "_")
my_comps3 <- list(c("T0_Yes", "T0_No"), c("T2_Yes", "T2_No"))


clono_plot<-ggplot(clono_vol_t_df_resp, aes(x= Timepoint_resp, y= Volume)) +
  geom_boxplot(aes(fill= Timepoint_resp), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Timepoint_resp), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("T0_Yes" = "#f58f71",
                              "T2_Yes" = "#f58f71",
                              "T0_No" ="#02dbc7",
                              "T2_No" ="#02dbc7"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("Unique clones")+
  xlab("")+
  ggtitle("Volume")+
  theme(plot.title = element_text(hjust=0.5)) +
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps3, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot
ggsave(plot = clono_plot, "Figures/Repertoire/T cell/Paper/Clone_volume_unique_responder_Tcell_timepoint.png", dpi = 300,  height = 6,  width = 6)


clono_vol_t_df_resp


summary_stats_T <- clono_vol_t_df_resp %>%
  group_by(Timepoint_resp) %>%
  summarise(
    Mean = mean(Volume),
    Median = median(Volume),
    Q1 = quantile(Volume, 0.25),
    Q3 = quantile(Volume, 0.75)
  )

summary_stats_T

####  diversity #######
####### Diversity choa ####

div_chao_t <- repDiversity(tcell$data, "chao1")
div_chao_t
div_chao_t<- as.data.frame(div_chao_t)

div_chao_t

div_chao_t$Sample <- rownames(div_chao_t)
rownames(div_chao_t)<- NULL

div_chao_t<- merge(div_chao_t, meta)

div_chao_t

div_chao_t<- div_chao_t[!div_chao_t$Timepoint == "T1", ]
div_chao_t<- div_chao_t[!div_chao_t$Timepoint == "T3", ]
div_chao_t$Status <- factor(div_chao_t$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0', 'Placebo T2'))


my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"),  c("PCV13 T0", "PCV13 T2"))


clono_plot<-ggplot(div_chao_t, aes(x= Status, y= Estimator)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  #ylab("Clones")+
  xlab("")+
  ggtitle("Chao 1")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_clones, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot

#c("Placebo T0", "Placebo T2"),
ggsave(plot = clono_plot, "Figures/Repertoire/T cell/Paper/Chao_diversity.png", dpi = 300,  height = 6,  width = 6)




####### Diversity simp ####

div_simp_t <- repDiversity(tcell$data, "inv.simp")
div_simp_t
div_simp_t<- as.data.frame(div_simp_t)

div_simp_t


meta
div_simp_t<- merge(div_simp_t, meta)

div_simp_t

div_simp_t<- div_simp_t[!div_simp_t$Timepoint == "T3", ]
div_simp_t<- div_simp_t[!div_simp_t$Timepoint == "T1", ]



div_simp_t$Status <- factor(div_simp_t$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0', 'Placebo T2'))


my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"), c("PCV13 T0", "PCV13 T2"))

colnames(div_simp_t)

clono_plot<-ggplot(div_simp_t, aes(x= Status, y= Value)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("1/D")+
  xlab("")+
  ggtitle("Inverse Simpson Index")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_clones, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot

#c("Placebo T0", "Placebo T2"),
ggsave(plot = clono_plot, "Figures/Repertoire/T cell/Paper/inver_simp_diversity.png", dpi = 300,  height = 6,  width = 6)

####### responder diversity

div_simp_t


div_simp_t_resp<- div_simp_t[!is.na(div_simp_t$Responder_final),]


div_simp_t_resp_T2<- div_simp_t_resp[div_simp_t_resp$Timepoint == "T2",]
table(div_simp_t_resp$Responder_final)


div_simp_t_resp$Timepoint_resp<- paste(div_simp_t_resp$Timepoint, div_simp_t_resp$Responder_final, sep = "_")
my_comps3 <- list(c("T0_Yes", "T0_No"), c("T2_Yes", "T2_No"))


clono_plot<-ggplot(div_simp_t_resp, aes(x= Timepoint_resp, y= Value)) +
  geom_boxplot(aes(fill= Timepoint_resp), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Timepoint_resp), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("T0_Yes" = "#f58f71",
                              "T2_Yes" = "#f58f71",
                              "T0_No" ="#02dbc7",
                              "T2_No" ="#02dbc7"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("1/D")+
  xlab("")+
  ggtitle("Inverse Simpson Index")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps3, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps4, map_signif_level = TRUE, test = "wilcox.test")

clono_plot

#c("Placebo T0", "Placebo T2"),
#ggsave(plot = clono_plot, "Figures/Repertoire/B cell/inver_simp_diversity.png", dpi = 300,  height = 6,  width = 6)
ggsave(plot = clono_plot, "Figures/Repertoire/T cell/Paper/inver_simp_diversity_overall_Responder_Tcell_timepoint.png", dpi = 300,  height = 6,  width = 6)

summary_stats <- div_simp_t_resp %>%
  group_by(Timepoint_resp) %>%
  summarise(
    Mean = mean(Value),
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

summary_stats
##########################
###### Rep clonality ####
##########################

# Set "homeo" to analyse relative abundance (also known as clonal space homeostasis),
#which is defined as the proportion of repertoire occupied by clonal groups with specific abundances..###

imm_hom_t <- repClonality(tcell$data, .method = "homeo")

imm_hom_t <- as.data.frame(imm_hom_t)
imm_hom_t

imm_hom_t$Sample <- rownames(imm_hom_t)
rownames(imm_hom_t)<- NULL

imm_hom_t
imm_hom_df_t<- merge(meta, imm_hom_t)

imm_hom_df_t


library(tidyr)

imm_hom_df_piv_t<-imm_hom_df_t %>% pivot_longer(cols=c('Rare (0 < X <= 1e-05)', 'Small (1e-05 < X <= 1e-04)', 'Medium (1e-04 < X <= 0.001)', 'Large (0.001 < X <= 0.01)', 'Hyperexpanded (0.01 < X <= 1)'),
                                                names_to='Clone_size',
                                                values_to='Values')


imm_hom_df_piv_t
imm_hom_df_piv_t<- imm_hom_df_piv_t[!imm_hom_df_piv_t$Timepoint == "T1", ]
imm_hom_df_piv_t<- imm_hom_df_piv_t[!imm_hom_df_piv_t$Timepoint == "T3", ]


imm_hom_df_piv_t$Status <- factor(imm_hom_df_piv_t$Status, levels=c('PCV13 T0', 'PCV13 T2', 'Placebo T0',  'Placebo T2'))


my_comps2 <- list(c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("Placebo T1", "Placebo T2"), c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"))
my_comps_clones <- list( c("Placebo T0", "Placebo T2"), c("PCV13 T0", "PCV13 T2"))

colnames(imm_hom_df_piv)



imm_hom_df_piv_t$Clone_size <- factor(imm_hom_df_piv_t$Clone_size, levels=c('Rare (0 < X <= 1e-05)', 'Small (1e-05 < X <= 1e-04)', 'Medium (1e-04 < X <= 0.001)', 'Large (0.001 < X <= 0.01)', 'Hyperexpanded (0.01 < X <= 1)'))


mycomps_1<- list( c("T0 PCV13", "T1 PCV13"), c("T0 PCV13", "T2 PCV13"), c("T0 PCV13", "T3 PCV13"), c("T0 Placebo", "T1 Placebo"), c("T0 Placebo", "T2 Placebo"), c("T0 Placebo", "T3 Placebo"))
my_comps_clones <-  list(c("PCV13 T0", "PCV13 T2"))
my_comps_clones_test<- list(c("T0", "T2"))


imm_hom_df_piv_t$Status <- as.factor(imm_hom_df_piv_t$Status)
imm_hom_df_piv_t$Treatment <- factor(imm_hom_df_piv_t$Treatment, levels=c('PCV13','Placebo'))




##### normality test 


# Split data by species
iris_split <- split(imm_hom_df_t$`Large (0.001 < X <= 0.01)`, imm_hom_df_t$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if greater than 0.05 then it is not normal 




my_comps_clones

#imm_hom_df_piv_t$Clone_size_status

my_comps_clones <- list( c("Rare (0 < X <= 1e-05)_Placebo_T0", "Rare (0 < X <= 1e-05)_Placebo_T2"), c("Rare (0 < X <= 1e-05)_PCV13_T0", "Rare (0 < X <= 1e-05)_PCV13_T2"))


class(imm_hom_df_piv_t$Clone_size_status)
imm_hom_df_piv_t$Clone_size_status<- as.factor(imm_hom_df_piv_t$Clone_size_status)


imm_hom_df_piv_t$Clone_size <- factor(imm_hom_df_piv_t$Clone_size, levels=c('Rare (0 < X <= 1e-05)', 'Small (1e-05 < X <= 1e-04)', 'Medium (1e-04 < X <= 0.001)', 'Large (0.001 < X <= 0.01)', 'Hyperexpanded (0.01 < X <= 1)'))

####### with facet ####### 
clono_plot<- ggplot(imm_hom_df_piv_t, aes(x= Clone_size, y= Values, fill=Status))+ 
  geom_boxplot(aes(fill= Status), outlier.shape = NA, shape =21, alpha = 0.4)+
  #geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 1) +
  geom_point(aes(fill=Status), shape= 21, position = position_jitterdodge(), size = 3)+
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("Percentage")+
  xlab("")+
  ggtitle("Clonal distribution")+
  theme(plot.title = element_text(hjust=0.5)) +
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  #geom_vline(xintercept = 3.5, linetype= "dotted", color = "black") +
  facet_wrap(vars(Treatment))+
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps2, method= "t.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps2, map_signif_level = TRUE, test = "wilcox.test") 

clono_plot
w <-clono_plot + geom_pwc(aes(group = Status), method = "t_test", label = "p.signif", hide.ns = TRUE)
w

ggsave("Clonality_homeo_stats_facet_t.png", w, path ="Figures/Repertoire/T cell/Paper/",  width = 20, height= 15, dpi = 300 )


################################################
################################################
#################### Top 10 ####################
################################################
################################################
library(immunarch)
library(tidyverse)
library(ggplot2)
library(circlize)
library(ggpubr)
#install.packages("wesanderson")
library(wesanderson)
#install.packages("ggalluvial")
library(ggalluvial)
library(dplyr)

#install.packages("viridis")  # Install
library(viridis) 
#install.packages("ggnewscale")
library(ggnewscale)


#################### B cell ####################



##### top 10 clones  ###
bcell_top10_list<-top(bcell$data, .n=100000000)
bcell_top10_list


top10_call<-do.call(rbind.data.frame, bcell_top10_list)

top10_call

top10_call$Sample <- rownames(top10_call)

top10_call

top10_call$Sample<- gsub("T0.*", "T0", top10_call$Sample)

top10_call$Sample<- gsub("T1.*", "T1", top10_call$Sample)

top10_call$Sample<- gsub("T2.*", "T2", top10_call$Sample)

top10_call$Sample<- gsub("T3.*", "T3", top10_call$Sample)

top10_call
meta
top10_call_meta<- merge(meta, top10_call, by="Sample")
top10_call_meta


top10_call_meta


### remove out of frame 
top10_call_meta<- top10_call_meta %>% filter(!CDR3.aa == "out_of_frame")
top10_call_meta

duplicates <- top10_call_meta[duplicated(top10_call_meta$CDR3.aa) | duplicated(top10_call_meta$CDR3.aa, fromLast = TRUE), ]
duplicates


write.csv(top10_call_meta, "Repertoire_csv/top10_call_meta_bcell.csv")

top10_call_meta<- read.csv("Repertoire_csv/top10_call_meta_bcell.csv")
top10_call_meta<-top10_call_meta[,-1]
top10_call_meta




####### split by vaccine and timepoint
### get the baseline top25 in both conditions together 


######

colnames(top10_call_meta)


PCV13<- top10_call_meta[top10_call_meta$Treatment %in% c( "PCV13") , ]

PCV13

PCV13_T0<- PCV13[PCV13$Timepoint %in% c( "T0") , ]

PCV13_T0


PCV13_T2<- PCV13[PCV13$Timepoint %in% c( "T2") , ]

PCV13_T2



Placebo<- top10_call_meta[top10_call_meta$Treatment %in% c( "Placebo") , ]

Placebo




Placebo_T0<- Placebo[Placebo$Timepoint %in% c( "T0") , ]

Placebo_T0

Placebo_T2<- Placebo[Placebo$Timepoint %in% c( "T2") , ]
Placebo_T2

### make a list of top10 start with PCV13 and Placebo 
library(dplyr)
colnames(PCV13)

PCV13_top10_clones<-PCV13 %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

PCV13_top10_clones



PCV13_top10_prop<-PCV13 %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

PCV13_top10_prop

###
Placebo

Placebo_count_2<- Placebo %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Placebo_count_2

duplicates <- Placebo_count_2[duplicated(Placebo_count_2$CDR3.aa) | duplicated(Placebo_count_2$CDR3.aa, fromLast = TRUE), ]
duplicates

#write.csv(Placebo_count_2, "Placebo_count2.csv")


#################
Placebo_count_2

Placebo_top10_Clones<-Placebo_count_2 %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

Placebo_top10_Clones


### for PCV

PCV13_count<- PCV13 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
PCV13_count

PCV13_count$CDR3.aa
# Find duplicated values in CDR3.aa
duplicates <- PCV13_count[duplicated(PCV13_count$CDR3.aa) | duplicated(PCV13_count$CDR3.aa, fromLast = TRUE), ]
duplicates
# Print duplicates



#write.csv(PCV13_count, "PCV_count.csv")


PCV13_top10_Clones<-PCV13_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

PCV13_top10_Clones

####PCV13 T0

PCV13_T0_count<- PCV13_T0 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
PCV13_T0_count

duplicates <- PCV13_T0_count[duplicated(PCV13_T0_count$CDR3.aa) | duplicated(PCV13_T0_count$CDR3.aa, fromLast = TRUE), ]
duplicates

#write.csv(PCV13_T0_count, "PCV_T0_count.csv")
#PCV13_T0_count<- read.csv( "PCV_T0_count.csv")


PCV13_T0_Clones<-PCV13_T0_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

PCV13_T0_Clones



PCV13_T0_prop<-PCV13_T0_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

PCV13_T0_prop

####PCV13 T2

PCV13_T2_count<- PCV13_T2 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
PCV13_T2_count

#write.csv(PCV13_T2_count, "PCV_T2_count.csv")
#PCV13_T2_count<- read.csv( "PCV_T2_count.csv")


PCV13_T2_Clones<-PCV13_T2_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
PCV13_T2_Clones

PCV13_T2_prop<-PCV13_T2_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

PCV13_T2_prop


PCV13_T2_prop25<-PCV13_T2_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:25) 

PCV13_T2_prop25



####Placebo T0

Placebo_T0_count<- Placebo_T0 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Placebo_T0_count

#write.csv(Placebo_T0_count, "Placebo_T0_count.csv")
#Placebo_T0_count<- read.csv("Placebo_T0_count.csv")



Placebo_T0_Clones<-Placebo_T0_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

Placebo_T0_Clones



Placebo_T0_prop<-Placebo_T0_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

Placebo_T0_prop


####Placebo T2

Placebo_T2_count<- Placebo_T2 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Placebo_T2_count

#write.csv(Placebo_T2_count, "Placebo_T2_count.csv")
#Placebo_T2_count<- read.csv("Placebo_T2_count.csv")


Placebo_T2_Clones<-Placebo_T2_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

Placebo_T2_Clones

Placebo_T2_prop<-Placebo_T2_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

Placebo_T2_prop


Placebo_T2_prop_25<-Placebo_T2_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:25) 

Placebo_T2_prop_25



#####################

######### take top 10 from T2 and follow them back to baseline  ###



PCV13_T2_count ### this is the all clones added together 
PCV13_T2_prop ### this is the top 10 

PCV13_T0_count


PCV13_T0_follow<- PCV13_T0_count[PCV13_T0_count$CDR3.aa %in% c("CMQALQTPHTF", "CQQSSKTPLTF", "CQQYSSLPYTF", "CASTVTGNNYYYGLDVW", "CHQYYDSPFTF", 
                                                               "CSSWDDSLSGRVF", "CFTRGGERGYSYGVYW", "CSSYTSSSSWVF", "CMQALQTPYTF", "CQAWDSSTVVF"),]


PCV13_T0_follow


PCV13_T2_prop$Timepoint<- "T2"
PCV13_T0_follow$Timepoint<- "T0"

colnames(PCV13_T2_prop)
colnames(PCV13_T0_follow)


PCV_follow<- rbind(PCV13_T2_prop, PCV13_T0_follow)


PCV_follow


########## try and add labels 

alluv <- ggplot(PCV_follow,
                aes(y = Proportion, x = Timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in PCV13") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion")  +
  theme(legend.text = element_text(size = 12))



alluv

ggsave("Alluv_PCV13_top10_T2_scale_withnumbers.png", alluv, path ="Figures/Repertoire/B cell/",   width = 15, height= 15, dpi = 300 )

alluv_p<- alluv + scale_fill_viridis_d(option =  "magma")

alluv_p

ggsave("Alluv_PCV13_top10_T2_scale_withnumbers_magma.png", alluv_p, path ="Figures/Repertoire/B cell/",   width = 15, height= 15, dpi = 300 )



###### Placebo ### 

Placebo_T2_prop

Placebo_T2_prop ### this is the top 10 

Placebo_T0_count


Placebo_T0_follow<- Placebo_T0_count[Placebo_T0_count$CDR3.aa %in% c("CQQSSKTPLTF", "CMQATQFPRTF", "CQQYNNWPPWTF", "CLQHNSYPWTF", "CQQYGSSPLTF", 
                                                                     "CASTVTGNNYYYGLDVW", "CQQYYSTPLTF", "CQEYNNWPPWTF", "CQAWDSSTVVF", "CQQSYSTPRTF"),]


Placebo_T0_follow


Placebo_T2_prop$Timepoint<- "T2"
Placebo_T0_follow$Timepoint<- "T0"

colnames(Placebo_T2_prop)
colnames(Placebo_T0_follow)


Placbo_follow<- rbind(Placebo_T2_prop, Placebo_T0_follow)


Placbo_follow


alluv <- ggplot(Placbo_follow,
                aes(y = Proportion, x = Timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in Placebo")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))



alluv


ggsave("Alluv_Placebo_top10_T2_scale_withnumbers.png", alluv, path ="Figures/Repertoire/B cell/",   width = 15, height= 15, dpi = 300 )


alluv_p<- alluv + scale_fill_viridis_d(option =  "mako")
alluv_p

ggsave("Alluv_Placebo_top10_T2_scale_withnumbers_mako.png", alluv_p, path ="Figures/Repertoire/B cell/",   width = 15, height= 15, dpi = 300 )


############### ############### ############### ############### 
############### T cells.       ############### ############### 
############### ############### ############### ############### 


##### top 10 clones  ###
tcell_top10_list<-top(tcell$data, .n=100000000)
tcell_top10_list


t_top10_call<-do.call(rbind.data.frame, tcell_top10_list)

t_top10_call

t_top10_call$Sample <- rownames(t_top10_call)

t_top10_call

t_top10_call$Sample<- gsub("T0.*", "T0", t_top10_call$Sample)

t_top10_call$Sample<- gsub("T1.*", "T1", t_top10_call$Sample)

t_top10_call$Sample<- gsub("T2.*", "T2", t_top10_call$Sample)

t_top10_call$Sample<- gsub("T3.*", "T3", t_top10_call$Sample)


t_top10_call_meta<- merge(meta, t_top10_call)



t_top10_call_meta

### remove out of frame 
t_top10_call_meta<- t_top10_call_meta %>% filter(!CDR3.aa == "out_of_frame")

write.csv(t_top10_call_meta, "Repertoire_csv/T_top10_call_meta_Tcell.csv")

t_top10_call_meta<- read.csv("Repertoire_csv/T_top10_call_meta_Tcell.csv")


####### Get the top counts for placebo and vaccine overall plus timepoints 


T_Placebo<- t_top10_call_meta[t_top10_call_meta$Treatment %in% c( "Placebo") , ]

T_Placebo

#write.csv(T_Placebo, "T_Placebo_count.csv")

T_Placebo_T0<- T_Placebo[T_Placebo$Timepoint %in% c( "T0") , ]

T_Placebo_T0



T_Placebo_T2<- T_Placebo[T_Placebo$Timepoint %in% c( "T2") , ]
T_Placebo_T2



###### PCV


T_PCV13<- t_top10_call_meta[t_top10_call_meta$Treatment %in% c( "PCV13") , ]

T_PCV13

T_PCV13_T0<- T_PCV13[T_PCV13$Timepoint %in% c( "T0") , ]

T_PCV13_T0


T_PCV13_T2<- T_PCV13[T_PCV13$Timepoint %in% c( "T2") , ]

T_PCV13_T2



##### Sum them up so that each sample is added to the total 

T_Placebo_count<- T_Placebo %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Placebo_count

T_placebo_top10_Clones<- T_Placebo_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_placebo_top10_Clones

T_placebo_top10_props<-T_Placebo_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 


T_placebo_top10_props


## PCV

T_PCV13_count<- T_PCV13 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_PCV13_count

T_PCV13_top10_Clones<-T_PCV13_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_PCV13_top10_Clones

T_PCV13_top10_props<-T_PCV13_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

T_PCV13_top10_props

### then get the summed totals for each timepoint ###


T_PCV13_T0_count<- T_PCV13_T0 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_PCV13_T0_count

#write.csv(T_PCV13_T0_count, "T_PCV_T0_count.csv")


T_PCV13_T0_Clones<-T_PCV13_T0_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_PCV13_T0_Clones



T_PCV13_T0_prop<-T_PCV13_T0_Clones %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
T_PCV13_T0_prop


####PCV13 T2

T_PCV13_T2_count<- T_PCV13_T2 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_PCV13_T2_count

#write.csv(T_PCV13_T2_count, "T_PCV_T2_count.csv")
#T_PCV13_T2_count<- read.csv("T_PCV_T2_count.csv")

T_PCV13_T2_Clones<-T_PCV13_T2_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_PCV13_T2_Clones

T_PCV13_T2_prop<-T_PCV13_T2_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

T_PCV13_T2_prop



####Placebo T0

T_Placebo_T0_count<- T_Placebo_T0 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Placebo_T0_count

#write.csv(T_Placebo_T0_count, "T_Placebo_T0_count.csv")


T_Placebo_T0_Clones<-T_Placebo_T0_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_Placebo_T0_Clones



T_Placebo_T0_prop<-T_Placebo_T0_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

T_Placebo_T0_prop


####Placebo T2

T_Placebo_T2_count<- T_Placebo_T2 %>% group_by(CDR3.aa) %>% summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Placebo_T2_count

#write.csv(T_Placebo_T2_count, "T_Placebo_T2_count.csv")


T_Placebo_T2_Clones<-T_Placebo_T2_count %>%
  arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_Placebo_T2_Clones

T_Placebo_T2_prop<-T_Placebo_T2_count %>%
  arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 

T_Placebo_T2_prop

######################################################################
######### take top 10 from T2 and follow them back to baseline  ###
######################################################################



T_PCV13_T0_count ### this is the all clones added together 
T_PCV13_T2_prop ### this is the top 10 


T_PCV13_T2_prop

T_PCV13_T0_follow<- T_PCV13_T0_count[T_PCV13_T0_count$CDR3.aa %in% c("CASTPGTDINQPQHF", "CSVEYRGAHEQYF", "C_ANHSVSSGSARQLTF", "CASSLGTDTQYF", "CAVGGYNFNKFYF", 
                                                                     "CASRTGTSDHEQYF", "CAISESSGWSQETQYF", "CASSQSQRGGYTF", "CAPPRARLMF", "CASSSTSTNTYEQYF"),]

T_PCV13_T0_follow


T_PCV13_T2_prop$Timepoint<- "T2"
T_PCV13_T0_follow$Timepoint<- "T0"

colnames(T_PCV13_T0_follow)
#T_PCV13_T2_prop<- T_PCV13_T2_prop[,-1]
colnames(T_PCV13_T2_prop)


T_PCV13_follow<- rbind(T_PCV13_T2_prop, T_PCV13_T0_follow)


T_PCV13_follow


alluv <- ggplot(T_PCV13_follow,
                aes(y = Proportion, x = Timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
           size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in PCV13") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))



alluv

ggsave("Alluv_PCV13_top10_T2_scale_withnumbers_Tcell.png", alluv, path ="Figures/Repertoire/T cell/",   width = 15, height= 15, dpi = 300 )

alluv_p<- alluv + scale_fill_viridis_d(option =  "magma")

alluv_p
ggsave("Alluv_PCV13_top10_T2_scale_withnumbers_magma_Tcell.png", alluv_p, path ="Figures/Repertoire/T cell/",   width = 15, height= 15, dpi = 300 )


###### Placebo ### 

T_Placebo_T0_count

T_Placebo_T2_prop ### this is the top 10 


T_Placebo_T0_follow<- T_Placebo_T0_count[T_Placebo_T0_count$CDR3.aa %in% c("CSVVWAKDTQYF", "CASSPTGGVGTEAFF", "CASSGADGGELFF", "CASSQEPHTDTQYF", "C_ANHSVSSGSARQLTF", 
                                                                           "CASSEFDSEKNNEQFF", "CASSLGTDTQYF", "CASSQERLARYGYTF", "CAGDSNYQLIW", "CAPSYSGGYNKLIF"),]


T_Placebo_T0_follow


T_Placebo_T2_prop$Timepoint<- "T2"
T_Placebo_T0_follow$Timepoint<- "T0"

colnames(T_Placebo_T2_prop)
colnames(T_Placebo_T0_count)


T_Placbo_follow<- rbind(T_Placebo_T2_prop, T_Placebo_T0_follow)


T_Placbo_follow



alluv <- ggplot(T_Placbo_follow,
                aes(y = Proportion, x = Timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in Placebo") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))


alluv

ggsave("Alluv_Placebo_top10_T2_scale_withnumbers_Tcell.png", alluv, path ="Figures/Repertoire/T cell/",   width = 15, height= 15, dpi = 300 )


alluv_p<- alluv + scale_fill_viridis_d(option =  "mako")
alluv_p
ggsave("Alluv_Placebo_top10_T2_scale_withnumbers_mako_Tcell.png", alluv_p, path ="Figures/Repertoire/T cell/",   width = 15, height= 15, dpi = 300 )


##############

################### responder ###############
###### B cell ########
top10_call_meta

#top10_call_meta$Responder_timepoint<- paste(top10_call_meta$Timepoint, top10_call_meta$Responder, sep = "_")
top10_call_meta$Responder_timepoint<- paste(top10_call_meta$Timepoint, top10_call_meta$Responder_final, sep = "_")

###### do this for the responders 

colnames(top10_call_meta)


Responder_Y<- top10_call_meta[top10_call_meta$Responder_final %in% c( "Yes") , ]

Responder_Y


Responder_N<- top10_call_meta[top10_call_meta$Responder_final %in% c( "No") , ]

Responder_N
###### and for timepoints so T0_Yes

Responder_Y_T0<- top10_call_meta[top10_call_meta$Responder_timepoint %in% c( "T0_Yes") , ]

Responder_Y_T0



Responder_Y_T2<- top10_call_meta[top10_call_meta$Responder_timepoint %in% c( "T2_Yes") , ]

Responder_Y_T2


### responder No timepoints
Responder_N_T0<- top10_call_meta[top10_call_meta$Responder_timepoint %in% c( "T0_No") , ]

Responder_N_T0




Responder_N_T2<- top10_call_meta[top10_call_meta$Responder_timepoint %in% c( "T2_No") , ]

Responder_N_T2



############## now do counts for the responders 

### first do total Responers and  non


Responder_Y_count<- Responder_Y %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Responder_Y_count

write.csv(Responder_Y_count, "Repertoire_csv/Responder_Y_count_update.csv")


Responder_Y_count

Responder_Y_top10_Clones<-Responder_Y_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

Responder_Y_top10_Clones

Responder_N_count<- Responder_N %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Responder_N_count

write.csv(Responder_N_count, "Repertoire_csv/Responder_N_count_update.csv")


Responder_N_count

Responder_N_top10_Clones<-Responder_N_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

Responder_N_top10_Clones

####


####Responder Y_T0

Resp_Y_T0_count<- Responder_Y_T0 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Resp_Y_T0_count

write.csv(Resp_Y_T0_count, "Repertoire_csv/Responder_Y_T0_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


Resp_Y_T0_Clones<-Resp_Y_T0_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
Resp_Y_T0_Clones



Resp_Y_T0_prop<-Resp_Y_T0_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
Resp_Y_T0_prop


#### T2
Resp_Y_T2_count<- Responder_Y_T2 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Resp_Y_T2_count

write.csv(Resp_Y_T2_count, "Repertoire_csv/Responder_Y_T2_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


Resp_Y_T2_Clones<-Resp_Y_T2_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
Resp_Y_T2_Clones



Resp_Y_T2_prop<-Resp_Y_T2_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
Resp_Y_T2_prop


#### responder NO
#### T0
Resp_N_T0_count<- Responder_N_T0 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Resp_N_T0_count

write.csv(Resp_N_T0_count, "Repertoire_csv/Responder_N_T0_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


Resp_N_T0_Clones<-Resp_N_T0_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
Resp_N_T0_Clones



Resp_N_T0_prop<-Resp_N_T0_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
Resp_N_T0_prop

#### T2
Resp_N_T2_count<- Responder_N_T2 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
Resp_N_T2_count

write.csv(Resp_N_T2_count, "Repertoire_csv/Responder_N_T2_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


Resp_N_T2_Clones<-Resp_N_T2_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
Resp_N_T2_Clones



Resp_N_T2_prop<-Resp_N_T2_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
Resp_N_T2_prop

#### bind them together try with proportion 

Resp_Y_T0_prop$Responder_timepoint<- "T0_Yes"
Resp_Y_T2_prop$Responder_timepoint<- "T2_Yes"



Responder_Y_timepoints_prop<- rbind(Resp_Y_T0_prop, Resp_Y_T1_prop, Resp_Y_T2_prop, Resp_Y_T3_prop)
Responder_Y_timepoints_prop

class(Responder_Y_timepoints_prop)


#########
####### want to find the top 10 clones in T2 responder and follow them back to T0 

Resp_Y_T2_prop ### look for the top10 here 

## in the T0 

Resp_Y_T0_count


Resp_Y_T0_follow<- Resp_Y_T0_count[Resp_Y_T0_count$CDR3.aa %in% c("CQQSSKTPLTF", "CASTVTGNNYYYGLDVW", "CHQYYDSPFTF", "CSSWDDSLSGRVF", "CMQALQTPYTF", 
                                                                  "CITKGRDITGQQVPRW", "CQQSYNSPRAF", "CATWDDSLRGVF", "CQVWDSRSNHVVF", "CQAWDSSTVVF"),]


Resp_Y_T0_follow

colnames(Resp_Y_T0_follow)
colnames(Resp_Y_T2_prop)

## add responder_timepoint column and bind them together 
#######
Resp_Y_T0_follow$Responder_timepoint<- "T0_Yes"

## bind them and use this for the alluvial 


Resp_yes_follow<- rbind(Resp_Y_T2_prop,Resp_Y_T0_follow )
Resp_yes_follow


Resp_yes_follow
alluv <- ggplot(Resp_yes_follow,
                aes(y = Proportion, x = Responder_timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in Responders") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))


alluv

alluv_p<- alluv + scale_fill_viridis_d(option =  "viridis")
alluv_p
ggsave("Alluv_top10_responder_Yes_update.png", alluv_p, path ="Figures/Repertoire/Responder/Bcell/",   width = 15, height= 15, dpi = 300 )



####### repeat for no responser ##

Resp_N_T2_prop


## in the T0 

Resp_N_T0_count


Resp_N_T0_follow<- Resp_N_T0_count[Resp_N_T0_count$CDR3.aa %in% c("CMQALQTPHTF", "CQQYSSLPYTF", "CQQSSKTPLTF", "CFTRGGERGYSYGVYW", "CAAWDDTLNGWVF", 
                                                                  "CQAWDSSTVVF", "CQQLNDYPLTF", "CMQALQTPLTF", "CAAWDDRLNGWVF", "CQQSYSTPRTF"),]


Resp_N_T0_follow

colnames(Resp_N_T0_follow)
colnames(Resp_N_T2_prop)

## add responder_timepoint column and bind them together 
#######
Resp_N_T2_prop$Responder_timepoint<- "T2_No"
Resp_N_T0_follow$Responder_timepoint<- "T0_No"

## bind them and use this for the alluvial 


Resp_no_follow<- rbind(Resp_N_T2_prop,Resp_N_T0_follow )


alluv <- ggplot(Resp_no_follow,
                aes(y = Proportion, x = Responder_timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in Non-Responders")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))


alluv

alluv_p<- alluv + scale_fill_viridis_d(option =  "mako")
alluv_p
ggsave("Alluv_top10_responder_No_update.png", alluv_p, path ="Figures/Repertoire/Responder/Bcell/",   width = 15, height= 15, dpi = 300 )


############################# T cells RESPONDER ####################

#### T cell ########
t_top10_call_meta

#top10_call_meta$Responder_timepoint<- paste(top10_call_meta$Timepoint, top10_call_meta$Responder, sep = "_")
t_top10_call_meta$Responder_timepoint<- paste(t_top10_call_meta$Timepoint, t_top10_call_meta$Responder_final, sep = "_")

###### do this for the responders 

colnames(top10_call_meta)
###### do this for the responders 



T_Responder_Y<- t_top10_call_meta[t_top10_call_meta$Responder_final %in% c( "Yes") , ]

T_Responder_Y


T_Responder_N<- t_top10_call_meta[t_top10_call_meta$Responder_final %in% c( "No") , ]

T_Responder_N

###### and for timepoints so T0_Yes

T_Responder_Y_T0<- t_top10_call_meta[t_top10_call_meta$Responder_timepoint %in% c( "T0_Yes") , ]

T_Responder_Y_T0

T_Responder_Y_T2<- t_top10_call_meta[t_top10_call_meta$Responder_timepoint %in% c( "T2_Yes") , ]

T_Responder_Y_T2

### responder No timepoints
T_Responder_N_T0<- t_top10_call_meta[t_top10_call_meta$Responder_timepoint %in% c( "T0_No") , ]

T_Responder_N_T0


T_Responder_N_T2<- t_top10_call_meta[t_top10_call_meta$Responder_timepoint %in% c( "T2_No") , ]

T_Responder_N_T2


############## now do counts for the responders 

### first do total Responers and  non


T_Responder_Y_count<- T_Responder_Y %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Responder_Y_count

write.csv(T_Responder_Y_count, "Repertoire_csv/Tcell/Responder_Y_count_update.csv")


T_Responder_Y_count

T_Responder_Y_top10_Clones<-T_Responder_Y_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_Responder_Y_top10_Clones

T_Responder_N_count<- T_Responder_N %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Responder_N_count

write.csv(T_Responder_N_count, "Repertoire_csv/Tcell/Responder_N_count_update.csv")


T_Responder_N_count

T_Responder_N_top10_Clones<-T_Responder_N_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 

T_Responder_N_top10_Clones

####


####Responder Y_T0

T_Resp_Y_T0_count<- T_Responder_Y_T0 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Resp_Y_T0_count

write.csv(T_Resp_Y_T0_count, "Repertoire_csv/Tcell/Responder_Y_T0_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


T_Resp_Y_T0_Clones<-T_Resp_Y_T0_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
T_Resp_Y_T0_Clones



T_Resp_Y_T0_prop<-T_Resp_Y_T0_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
T_Resp_Y_T0_prop


#### T2
T_Resp_Y_T2_count<- T_Responder_Y_T2 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Resp_Y_T2_count

write.csv(T_Resp_Y_T2_count, "Repertoire_csv/Tcell/Responder_Y_T2_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


T_Resp_Y_T2_Clones<-T_Resp_Y_T2_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
T_Resp_Y_T2_Clones



T_Resp_Y_T2_prop<-T_Resp_Y_T2_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
T_Resp_Y_T2_prop


#### responder NO
#### T0
T_Resp_N_T0_count<- T_Responder_N_T0 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Resp_N_T0_count

write.csv(T_Resp_N_T0_count, "Repertoire_csv/Tcell/Responder_N_T0_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


T_Resp_N_T0_Clones<-T_Resp_N_T0_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
T_Resp_N_T0_Clones



T_Resp_N_T0_prop<-T_Resp_N_T0_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
T_Resp_N_T0_prop


#### T2
T_Resp_N_T2_count<- T_Responder_N_T2 %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones= sum(Clones), Proportion= sum(Proportion))
T_Resp_N_T2_count

write.csv(T_Resp_N_T2_count, "Repertoire_csv/Tcell/Responder_N_T2_count_update.csv")
#Resp_Y_T0_count<- read.csv("Repertoire_csv/Responder_Y_T0_count.csv")


T_Resp_N_T2_Clones<-T_Resp_N_T2_count %>%
  dplyr::arrange(desc(Clones)) %>%
  dplyr::slice(1:10) 
T_Resp_N_T2_Clones



T_Resp_N_T2_prop<-T_Resp_N_T2_count %>%
  dplyr::arrange(desc(Proportion)) %>%
  dplyr::slice(1:10) 
T_Resp_N_T2_prop

##### want to find the top 10 clones in T2 responder and follow them back to T0 

T_Resp_N_T2_prop ### look for the top10 here 

## in the T0 

T_Resp_N_T0_count


T_Resp_N_T0_follow<- T_Resp_N_T0_count[T_Resp_N_T0_count$CDR3.aa %in% c("CASTPGTDINQPQHF", "CAVGGYNFNKFYF", "CASRTGTSDHEQYF", "CASSQSQRGGYTF", "CAPPRARLMF", 
                                                                        "C_ANHSVSSGSARQLTF", "CASSYLERTGQPQHF", "CASSPSNYGYTF", "CASSLSRPSGSPLHF", "CASSRDFGETQYF"),]


T_Resp_N_T0_follow

colnames(T_Resp_N_T0_follow)
colnames(T_Resp_N_T2_prop)

## add responder_timepoint column and bind them together 
#######
T_Resp_N_T0_follow$Responder_timepoint<- "T0_No"
T_Resp_N_T2_prop$Responder_timepoint<- "T2_No"

## bind them and use this for the alluvial 


T_Resp_no_follow<- rbind(T_Resp_N_T2_prop,T_Resp_N_T0_follow )

T_Resp_no_follow


T_Resp_no_follow


alluv <- ggplot(T_Resp_no_follow,
                aes(y = Proportion, x = Responder_timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in Non-Responders")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))

alluv

alluv_p<- alluv + scale_fill_viridis_d(option =  "mako")
alluv_p

ggsave("Alluv_top10_responder_No_update.png", alluv_p, path ="Figures/Repertoire/Responder/Tcell/",   width = 15, height= 15, dpi = 300 )


##########

###### want to find the top 10 clones in T2 responder and follow them back to T0 

T_Resp_Y_T2_prop ### look for the top10 here 

## in the T0 

T_Resp_Y_T0_count


T_Resp_Y_T0_follow<- T_Resp_Y_T0_count[T_Resp_Y_T0_count$CDR3.aa %in% c("CSVEYRGAHEQYF", "CAISESSGWSQETQYF", "CASSSTSTNTYEQYF", "CASSLGTDTQYF", "CASSLRVASPAYNEQFF", 
                                                                        "CASSLAAGANEQFF", "C_ANHSVSSGSARQLTF", "CASSPGTSRTPETQYF", "CAASPPESGGYNKLIF", "CSVEGTGGYEQYF"),]


T_Resp_Y_T0_follow

colnames(T_Resp_Y_T0_follow)
colnames(T_Resp_Y_T0_count)

## add responder_timepoint column and bind them together 
#######
T_Resp_Y_T0_follow$Responder_timepoint<- "T0_Yes"
T_Resp_Y_T2_prop$Responder_timepoint<- "T2_Yes"

## bind them and use this for the alluvial 


T_Resp_yes_follow<- rbind(T_Resp_Y_T2_prop,T_Resp_Y_T0_follow )


T_Resp_yes_follow

##### alluvial #####

Resp_yes_follow
alluv <- ggplot(T_Resp_yes_follow,
                aes(y = Proportion, x = Responder_timepoint, stratum = CDR3.aa, alluvium = CDR3.aa, fill = CDR3.aa)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray") +
  geom_stratum(aes(fill = CDR3.aa), alpha = 0.5) +
  geom_alluvium(aes(fill = CDR3.aa)) +
  geom_text(stat = "flow", aes(label = sprintf("%.2f", Proportion)), 
            size = 7, color = "black", fontface = "bold", check_overlap = TRUE, na.rm = TRUE) + 
  theme_bw() +
  theme(legend.position = "right") +
  ggtitle("Top 10 Clones in Responders") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Timepoint") +
  ylab("Proportion") +
  theme(legend.text = element_text(size = 12))


alluv

alluv_p<- alluv + scale_fill_viridis_d(option =  "viridis")
alluv_p

ggsave("Alluv_top10_responder_Yes_update.png", alluv_p, path ="Figures/Repertoire/Responder/Tcell/",   width = 15, height= 15, dpi = 300 )


