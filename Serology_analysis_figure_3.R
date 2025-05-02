###### Vaciriss Serology ######
## boxplots over cohort time ##
## Fold change: bubble plot and FC correlation  ##
## Volcano plots  ##


#####  responder at the bottom? ###


library(tidyverse)
library(ggplot2)
library(readr)
library(limma)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(viridis)
library(plyr)
library(dplyr)

###### load in data 

getwd()

######################## serology boxplots ################################## 


WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")
WHO_T0_T2
WHO_T0_T2<- WHO_T0_T2[,-c(1,2)]

WHO_T0_T2$Status<- factor(WHO_T0_T2$Status, levels = c( "PCV13 T0", "PCV13 T2", "Placebo T0", "Placebo T2")) 


my_comps<- list(c("Placebo T0", "Placebo T2"), c("PCV13 T0", "PCV13 T2"), c("Placebo T0", "PCV13 T0"), c("Placebo T2", "PCV13 T0") )
my_comps_2<- list( c("PCV13 T0", "PCV13 T2"))
my_comps_3<- list( c("PCV13 T0", "PCV13 T2"), c("Placebo T0", "Placebo T2"))

WHO_T0_T2

colnames(WHO_T0_T2)

table(WHO_T0_T2$Status)



# Split data by species
iris_split <- split(WHO_T0_T2$IgG_type.23F, WHO_T0_T2$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 


################### boxplot Serotpye with cohort_time ######

serotype<-ggplot(WHO_T0_T2, aes(x= Status, y= IgG_type.3)) +
  geom_boxplot(aes(fill= Status), width=0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              "Placebo T2" ="#577a6c",
                              "PCV13 T0" = "#f3c6d9",
                              "PCV13 T2" = "#d32c75"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("IgG ug/ml")+
  xlab("")+
  ggtitle("3")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") +
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = list( c("PCV13 T0", "PCV13 T2")), method= "annova",  label = "p.signif", hide.ns = TRUE)
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps, label = "p.signif", hide.ns = TRUE)
  geom_signif(comparisons = my_comps_2, map_signif_level = TRUE, test = "wilcox.test")

serotype

#c("Placebo T0", "Placebo T2"),
ggsave(plot = serotype, "Figures/Serology/Serotype23F_box.png", dpi = 300,  height = 6,  width = 6)
########
########

colnames(WHO_T0_T2)
WHO_T0_T2

summary <- WHO_T0_T2 %>%
  group_by(Status) %>%
  dplyr::summarise(
    median_percentage = median(IgG_type.23F, na.rm = TRUE),
    Q1 = quantile(IgG_type.23F, 0.25, na.rm = TRUE),
    Q3 = quantile(IgG_type.23F, 0.75, na.rm = TRUE),
    sample_count=n()
  )

summary




##### fold change 


WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")


WHO_T0_T2
WHO_T0_T2<- WHO_T0_T2[,-c(1,2)]

WHO_T0_T2$group

Serology_merge<-WHO_T0_T2
Serology_merge

PCV<- Serology_merge[Serology_merge$group== "PCV13", ]
PCV

nrow(PCV)

PCV_T0<- PCV[PCV$timepoint=="T0",]
PCV_T2<- PCV[PCV$timepoint=="T2",]

nrow(PCV_T0) #90 
nrow(PCV_T2) #48

Placebo<- Serology_merge[Serology_merge$group== "Placebo", ]


Placebo_T0<- Placebo[Placebo$timepoint=="T0",]
Placebo_T2<- Placebo[Placebo$timepoint=="T2",]

nrow(Placebo_T0) #93
nrow(Placebo_T2) #38

Serology_merge



# Replace NA in all IgG columns with 0.15 (lowest detectable value)
data <- Serology_merge %>%
  dplyr::mutate(across(starts_with("IgG_type"), ~ ifelse(is.na(.), 0.15, .)))

data


# Calculate the difference between timepoints for each IgG type and each patient
df_long <- data %>%
  pivot_longer(cols = starts_with("IgG_type"), names_to = "IgG_type", values_to = "IgG_value") %>%
  group_by(patient_num, IgG_type)

df_long

#df_long$oldID<- df_long$patient_num
#df_long

df_FC <- df_long %>%
  pivot_wider(
    id_cols = c(patient_num, IgG_type, group, Responder_final),
    names_from = timepoint,
    values_from = IgG_value
  ) %>%
  select(patient_num, IgG_type, T0, T2, group, Responder_final)

df_FC


df_FC$FC_T2<- df_FC$T2/df_FC$T0

df_FC


write.csv(df_FC, "WHO_Serology/Fold_change_df_anon2.csv")


################ bubble plot of fold change ################# 



df_FC <- read.csv("WHO_Serology/Fold_change_df_anon2.csv")
df_FC


df_FC

table(df_FC$IgG_type)

mean_fc_per_serotype <- df_FC %>%
  group_by(IgG_type, group) %>%
  dplyr::summarise(mean_FC_T2 = mean(FC_T2, na.rm = TRUE), .groups = "drop")
# View the result
mean_fc_per_serotype


mean_fc_per_serotype_df<- as.data.frame(mean_fc_per_serotype)

mean_fc_per_serotype_df

wide_format <- mean_fc_per_serotype_df %>%
  pivot_wider(
    names_from = group,
    values_from = mean_FC_T2,
    names_prefix = "FC_"
  )

wide_format



wide_format$Ratio<- wide_format$FC_PCV13/wide_format$FC_Placebo

wide_format

# Pivot the data longer
long_format <- wide_format %>%
  pivot_longer(
    cols = starts_with("FC_") | starts_with("Ratio"),
    names_to = "Metric",
    values_to = "Value"
  )

long_format$Metric

long_format$Metric<-gsub("FC_PCV13","PCV13",as.character(long_format$Metric))
long_format$Metric<-gsub("FC_Placebo","Placebo",as.character(long_format$Metric))

long_format$IgG_type<-gsub("IgG_type.4","4",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.23F","23F",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.6A","6A",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.7F","7F",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.5","5",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.14","14",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.1","1",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.6B","6B",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.19A","19A",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.18C","18C",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.19F","19F",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.9V","9V",as.character(long_format$IgG_type))
long_format$IgG_type<-gsub("IgG_type.3","3",as.character(long_format$IgG_type))

long_format

long_format_ratio<- long_format[long_format$Metric== "Ratio", ]


long_format$IgG_type<- factor(long_format$IgG_type, levels = c("14","18C", "6A", "23F", "4", "1", 
                                                               "19A", "6B", "5", "9V", "7F", "19F", "3"))


heat <- ggplot(long_format, aes(x = Metric, y = IgG_type, color= Value)) +
  geom_point(aes( size= Value), na.rm = FALSE)+ 
  labs(title = "Average IgG value",
       x = "Treatment",
       y = "Strain") +
  theme_bw() +
  scale_y_discrete(limits = rev(levels(factor(long_format$IgG_type))))+
  scale_size(range = c(0.1, 25)) 


map<- heat+ scale_color_viridis_c(option = "inferno")
map

ggsave(plot = map, "Figures/Serology/Bubble_ratio_PCV_placebo_new.png", dpi = 300, height = 12, width =12, units = "in")




########## correlation of FC changes ##########################
library(corrplot)
library(ggcorrplot)

df_FC <- read.csv("WHO_Serology/Fold_change_df_anon.csv")
df_FC

FC_matrix<- df_FC %>% select(patient_num , IgG_type, FC_T2, Treatment)


FC_matrix
FC_matrix_PCV<- FC_matrix[FC_matrix$Treatment== "PCV13", ]

FC_matrix_PCV<- na.omit(FC_matrix_PCV)

FC_matrix_pivot <- FC_matrix_PCV %>%
  pivot_wider(
    names_from = IgG_type,  # Make each Serotype a column
    values_from = FC_T2,       # Fill the values of those columns with the FC values
    values_fill = list(FC = NA)  # Fill NAs if needed
  )

FC_matrix_pivot

names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.1"] <- "1"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.3"] <- "3"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.4"] <- "4"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.5"] <- "5"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.6A"] <- "6A"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.6B"] <- "6B"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.7F"] <- "7F"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.9V"] <- "9V"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.14"] <- "14"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.18C"] <- "18C"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.19A"] <- "19A"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.19F"] <- "19F"
names(FC_matrix_pivot)[names(FC_matrix_pivot) == "IgG_type.23F"] <- "23F"

FC_matrix_pivot_2<- FC_matrix_pivot %>% column_to_rownames(var="oldID.x")

FC_matrix_pivot_2<- FC_matrix_pivot_2 %>% select(!Treatment  )

FC_matrix_pivot_2
FC_matrix_pivot_2<- as.matrix(FC_matrix_pivot_2)


FC_matrix_pivot_2


# Create the correlation matrix
cor_matrix <- cor(FC_matrix_pivot_2, use = "pairwise.complete.obs")
cor_matrix


# Plot the correlation matrix using ggcorrplot


corrplot(cor_matrix, type = 'lower', order = 'hclust', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45, col = rev(COL2('RdBu', 10)))


# Create a PDF file to save the plot
pdf("Figures/Serology/FC_corrplot_update.pdf", width = 10, height = 10)  # Specify the file name and dimensions in inches

corrplot(cor_matrix, type = 'lower', order = 'hclust', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45, col = rev(COL2('RdBu', 10)))

# Close the PDF device to save the plot
dev.off()



################ ################ ################ 
################ volcano plots ################# 
################ ################ ################ 

WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")

WHO_T0_T2<- WHO_T0_T2[,-1]

WHO_T0_T2


#replace NA with lower limit 0.15

WHO_T0_T2[is.na(WHO_T0_T2)] <- 0.15


WHO_T0_T2


WHO_T0_T2
#### filter samples based on timepoint 
df_norm_T0 <- filter(WHO_T0_T2, timepoint == "T0") 
df_norm_T0
df_norm_T2 <- filter(WHO_T0_T2, timepoint == "T2") 


#Log the data 

##### log T0 and T2 ###
df_log_T0 <- df_norm_T0 %>% mutate_at(vars(2:14), log2) %>% 
  select(2:14, anon_sample_ID, patient_num, timepoint, group, Status) %>% 
  as.data.frame()


df_log_T2 <- df_norm_T2 %>% mutate_at(vars(2:14), log2) %>% 
  select(2:14, anon_sample_ID, patient_num, timepoint, group, Status) %>% 
  as.data.frame()
df_log_T2


############################################################
###select data for comparison and remove rows with NA's ###

df_log_T0 <- df_log_T0 %>%  drop_na()

df_log_T2 <- df_log_T2 %>%  drop_na()




#select columns for analysis for limma to create volcano plots 
#WHO

analysis_T0 <- df_log_T0 %>% select("anon_sample_ID", 1:13)
analysis_T2 <- df_log_T2 %>% select("anon_sample_ID", 1:13)



#make matrix of values, samples as columns and frequencies as rows

###### T0 and T2
count_comparison_tmydf_T0 = setNames(data.frame(t(analysis_T0[,-1])), analysis_T0[,1])%>% as.matrix() 
str(count_comparison_tmydf_T0)
dim(count_comparison_tmydf_T0)


count_comparison_tmydf_T2 = setNames(data.frame(t(analysis_T2[,-1])), analysis_T2[,1])%>% as.matrix() 
str(count_comparison_tmydf_T2)
dim(count_comparison_tmydf_T2)


######### creat factor levels ####


f_T0 <- factor(df_log_T0$group, levels=c("Placebo", "PCV13"))
f_T2 <- factor(df_log_T2$group, levels=c("Placebo", "PCV13"))



design_T0 <- model.matrix(~0+f_T0)
design_T2 <- model.matrix(~0+f_T2)


colnames(design_T0) <- c("Placebo", "PCV13")
colnames(design_T2) <- c("Placebo", "PCV13")


#To make all pair-wise comparisons between the three groups one could proceed

fit_T0 <- lmFit(count_comparison_tmydf_T0, design_T0)
fit_T2 <- lmFit(count_comparison_tmydf_T2, design_T2)


##########


contrast.matrixT0 <- makeContrasts(PCV13-Placebo,
                                   levels=design_T0)
contrast.matrixT2 <- makeContrasts(PCV13-Placebo,
                                   levels=design_T2)


#####
###############

fit2_T0 <- contrasts.fit(fit_T0, contrast.matrixT0)
fit2_T0 <- eBayes(fit2_T0)
topTable(fit2_T0)

fit2_T2 <- contrasts.fit(fit_T2, contrast.matrixT2)
fit2_T2 <- eBayes(fit2_T2)
topTable(fit2_T2)


df_statsT0 <- topTable(fit2_T0, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("PNC")
df_statsT2 <- topTable(fit2_T2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("PNC")



write.csv(df_statsT0, "WHO_Serology/Limma_WHO_T0_PCVvPlacebo_anon.csv")


############### load in the limma files again and run the volcanos 
df<- read.csv("WHO_Serology/Limma_WHO_T0_PCVvPlacebo_anon.csv")




# add a column of NAs
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$logFC > 0.5 & df$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$logFC < -0.5 & df$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df$delabel <- NA
df$delabel<- df$PNC
##### below works to cut out the names of those that aren't significant
df$delabel[df$diffexpressed != "NO"] <- df$PNC[df$diffexpressed != "NO"]

mycolors <- c("#d102ce", "firebrick1", "grey")
names(mycolors) <- c("UP", "DOWN", "NO")

VP <- ggplot(data=df, aes(x=logFC, y=-log10( adj.P.Val), label=delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 4)+
  theme_classic()+
  ggtitle("T0 Placebo V PCV13")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold", size = 16))+
  geom_vline(xintercept=c(-0.5, 0.5),linetype = "dashed",  col="black") +
  geom_hline(yintercept=-log10(0.05),linetype = "dashed",  col="black")+
  scale_fill_manual(values=mycolors)+
  geom_label_repel(size = 3, force = 3, max.overlaps =25)+
  xlim(c(-2,2))+
  ylim(c(0,2.5))


VP

ggsave(plot = VP, "Figures/Serology/VP_WHO_T0_PCVvPlac_PCV.png", dpi = 300, height = 5, width =5, units = "in")


############# ############# ############# #############

############# Serology Responders #############




WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")
WHO_T0_T2
nrow(WHO_T0_T2)
WHO_T0_T2<- WHO_T0_T2[,-c(1,2)]

colnames(WHO_T0_T2)
Serology_merge<- WHO_T0_T2


PCV<- Serology_merge[Serology_merge$group== "PCV13", ]
PCV


# Replace NA in all IgG columns with 0.15
data <- PCV %>%
  dplyr::mutate(across(starts_with("IgG_type"), ~ ifelse(is.na(.), 0.15, .)))

data


# Calculate the difference between timepoints for each IgG type and each patient
df_long <- PCV %>%
  pivot_longer(cols = starts_with("IgG_type"), names_to = "IgG_type", values_to = "IgG_value") %>%
  group_by(patient_num, IgG_type)

df_long

df_T0<- df_long[df_long$timepoint == "T0",]

df_T0
nrow(df_T0)

unique(df_T0$patient_num)

colnames(df_T0)[which(names(df_T0) == "IgG_value")] <- "Baseline"

colnames(df_T0)
df_T0<- df_T0 %>% select(patient_num,  IgG_type, Baseline)
df_T0



df_T2<- df_long[df_long$timepoint == "T2",]

df_T2
unique(df_T2$patient_num)


colnames(df_T2)[which(names(df_T2) == "IgG_value")] <- "T2"
colnames(df_T0)

df_T2<- df_T2 %>% select(patient_num,  IgG_type, T2)
df_T2

#df_FC<- merge(df_T0, df_T1, by = c("PIN", "patient_ID", "IgG_type"), all = TRUE)
#df_FC

df_FC<- merge(df_T0, df_T2, by = c("patient_num",  "IgG_type"), all = TRUE)

df_FC
colnames(df_FC)
df_FC



#df_FC$FC_T1 <- df_FC$T1/df_FC$Baseline

df_FC$FC_T2<- df_FC$T2/df_FC$Baseline


df_FC
write.csv(df_FC, "WHO_Serology/Fold_change_df_anon_2.csv")

df_FC <- read.csv("WHO_Serology/Fold_change_df_anon_2.csv")
df_FC



##### repeat for serotype 23F



Ig23F<- df_FC[df_FC$IgG_type == "IgG_type.23F",]


Ig23F


###### find the average of the IgG4 for T1 and T2 

median(Ig23F$Baseline, na.rm = TRUE) #  

mean(Ig23F$T2, na.rm = TRUE) # 
median(Ig23F$T2, na.rm = TRUE) #  


Ig23F

#####T2 rule FC change of >2 & Post value is more than the average of baseline it is a responder,  or if no FC is present post value is 2x the average value of T2


Ig23F$Responder_T2<- "No"
Ig23F

Ig23F$Responder_T2[Ig23F$FC_T2 > 2 ] <- "Yes" #& Ig23F$T2 > 0.7

#Ig23F$Responder_T2[Ig23F$Responder_T2== "No" & Ig23F$T2 > 12.2 &  Ig23F$T2 > Ig23F$Baseline] <- "Yes"

Ig23F

table(Ig23F$Responder_T2)


##### repeat for serotype 6A



Ig6A<- df_FC[df_FC$IgG_type == "IgG_type.6A",]


Ig6A

Ig6A$Responder_T2<- "No"
Ig6A



Ig6A$Responder_T2[Ig6A$FC_T2 > 2 ] <- "Yes" # & Ig6A$T2 > mean_baseline

#Ig6A$Responder_T2[Ig6A$Responder_T2== "No" & Ig6A$T2 > 15.8 &  Ig6A$T2 > Ig6A$Baseline] <- "Yes"

Ig6A

table(Ig6A$Responder_T2)



##### repeat for serotype 7F



Ig7F<- df_FC[df_FC$IgG_type == "IgG_type.7F",]


Ig7F



Ig7F$Responder_T2<- "No"
Ig7F

Ig7F$Responder_T2[Ig7F$FC_T2 > 2 ] <- "Yes" #& Ig7F$T2 > mean_baseline


Ig7F

table(Ig7F$Responder_T2)

#####repeat for serotype 14



Ig14<- df_FC[df_FC$IgG_type == "IgG_type.14",]


Ig14



Ig14$Responder_T2<- "No"
Ig14



Ig14$Responder_T2[Ig14$FC_T2 > 2 ] <- "Yes" # & Ig6A$T2 > mean_baseline


Ig14

table(Ig14$Responder_T2)

##### repeat for serotype 18C



Ig18C<- df_FC[df_FC$IgG_type == "IgG_type.18C",]


Ig18C



Ig18C$Responder_T2<- "No"
Ig18C



Ig18C$Responder_T2[Ig14$FC_T2 > 2 ] <- "Yes" # & Ig6A$T2 > mean_baseline


Ig18C

table(Ig18C$Responder_T2)


#### rbind them



resp<- rbind(Ig14, Ig18C, Ig6A ,Ig23F)

resp


# Create a new dataframe with the desired structure
library(dplyr)
resp

##### if 


df_responder <- resp %>%
  group_by(patient_num) %>%
  dplyr::summarize(
    Responder_T2 = ifelse(sum(Responder_T2 == "Yes") >= 2, "Yes", "No"),
    Baseline_present = ifelse(any(is.na(Baseline)), "No", "Yes"),
    T2_present = ifelse(any(is.na(T2)), "No", "Yes")
  ) %>%
  ungroup()

df_responder

df_responder <- resp %>%
  group_by(patient_num) %>%
  dplyr::summarize(
    Responder_T2 = ifelse(sum(Responder_T2 == "Yes") >= 2, "Yes", "No"),
    Baseline_present = ifelse(any(is.na(Baseline)), "No", "Yes"),
    T2_present = ifelse(any(is.na(T2)), "No", "Yes")
  ) %>%
  ungroup() %>%
  mutate(
    T2_adjusted = ifelse(T2_present == "No", NA, Responder_T2)
  )

df_responder

# View the result
print(df_responder)

df_responder
# View the result
print(df_responder)

df_responder
df_responder$
# View the result
df_responder
table(df_responder$Responder_final)

write.csv(df_responder, "WHO_Serology/Responder_T2_anon_patientID.csv")
df_responder<- read.csv("WHO_Serology/Responder_T2_anon_patientID.csv")

table(df_responder$T2_adjusted)
df_responder

################### quick check 

#colnames(df_responder)[which(names(df_responder) == "Responder_final")] <- "Responder_final_test"
df_responder


WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")

WHO_T0_T2

colnames(WHO_T0_T2)
colnames(df_responder)

 
merger<- merge(WHO_T0_T2, df_responder, by= "patient_num")
 
colnames(merger)
test <- merger[ !is.na(merger$Responder_final) & 
                  !is.na(merger$T2_adjusted) & 
                  merger$Responder_final != merger$T2_adjusted, ]
test


##############################################################################
################## responder serology boxplots responder v non ###############



WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")

colnames(WHO_T0_T2)

serology_long_resp_filter <- WHO_T0_T2 %>%
  filter(!is.na(Responder_final))



serology_long_resp_filter$time_resp<- paste(serology_long_resp_filter$timepoint, serology_long_resp_filter$Responder_final, sep= "_")

serology_long_resp_filter$time_resp

serology_long_resp_filter$time_resp<- factor(serology_long_resp_filter$time_resp, levels = c("T0_No", "T2_No", "T0_Yes", "T2_Yes"))

colnames(serology_long_resp_filter)
# Pivot the data longer
serology_long_resp_filter <- serology_long_resp_filter %>%
  pivot_longer(
    cols = starts_with("IgG_type") ,
    names_to = "IgG_type",
    values_to = "IgG_value"
  )

serology_long_resp_filter

IgG23F_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.23F", ]
IgG19A_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.19A", ]
IgG5_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.5", ]
IgG14_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.14", ]
IgG1_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.1", ]
IgG4_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.4", ]

IgG19F_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.19F", ]
IgG18C_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.18C", ]
IgG7F_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.7F", ]
IgG9V_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.9V", ]
IgG6B_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.6B", ]
IgG6A_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.6A", ]
IgG3_resp<- serology_long_resp_filter[serology_long_resp_filter$IgG_type== "IgG_type.3", ]


#####

my_comps_1<- list(c("T0_No", "T2_No"), c("T0_Yes", "T2_Yes"))

# Split data by species
iris_split <- split(IgG3_resp$IgG_value, IgG3_resp$time_resp)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if greater than 0.05 then it is not normal 


colnames(IgG1_resp)
t_test_results <- compare_means(IgG_value ~ time_resp, data = IgG3_resp, method = "wilcox.test")
t_test_results

significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results


cyto_plot<-ggplot(IgG23F_resp, aes(x= time_resp, y= IgG_value)) +
  geom_boxplot(aes(fill= time_resp), width= 0.35, outlier.shape = NA, alpha = 0.4)+
  #geom_hline(yintercept = 54:68, linetype= "dashed", color = "grey",  alpha = 0.6)+
  #geom_jitter(aes(fill= Status, shape= Responder), width=0.25, alpha=1, size = 3) +
  geom_jitter(aes(fill= time_resp),  shape= 21, width=0.25, alpha=1, size = 3) +
  #scale_shape_manual(values = c("Yes" = 21, "No" = 24))+
  scale_fill_manual(values= c( "T0_No" = "#02dbc7", 
                               "T1_No" = "#02b5a4",
                               "T2_No" = "#009688",
                               "T0_Yes" = "#f58f71",
                               "T1_Yes" = "#f57953",
                               "T2_Yes" = "#FF7043"))+
  #scale_fill_manual(values = c("PS" = "green", 
  #                             "YG"= "blue"))+
  #scale_shape()+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black")
  ) +
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  #annotation_ticks(sides = "l", type = "minor", outside = TRUE)+
  ylab("IgG ug/ml")+
  xlab("")+
  ggtitle("23F")+
  #ylim(0, 100) +
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  #geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black")  +
  #geom_vline(xintercept = 1.5, linetype= "dotted", color = "black")  +
  #geom_hline(yintercept = 54:68, linetype= "dashed", color = "red")  +
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_1, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE)# + #label.y = c(95, 90, 85, 80, 75, 40) 
  #scale_y_continuous(trans = "log1p")
  geom_signif(comparisons = my_comps_1, map_signif_level = TRUE, test = "wilcox.test")


cyto_plot

ggsave(plot = cyto_plot, "Figures/Responder/Serology/23F_responder.png", dpi = 300,  height = 6,  width = 6)

