############### Vaciriss Figure 3 ##########
## serology  ##
####### Volcano
##### should add the fold change bubble plot here for the upload
## fold change graphs, heatmap and age/sex boxplots 
##  Fold change correlation plts 

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
  dplyr::select(2:14, anon_sample_ID, patient_num, timepoint, group, Status) %>% 
  as.data.frame()


df_log_T2 <- df_norm_T2 %>% mutate_at(vars(2:14), log2) %>% 
  dplyr::select(2:14, anon_sample_ID, patient_num, timepoint, group, Status) %>% 
  as.data.frame()
df_log_T2


############################################################
###select data for comparison and remove rows with NA's ###

df_log_T0 <- df_log_T0 %>%  drop_na()

df_log_T2 <- df_log_T2 %>%  drop_na()




#select columns for analysis for limma to create volcano plots 
#WHO

analysis_T0 <- df_log_T0 %>% dplyr::select("anon_sample_ID", 1:13)
analysis_T2 <- df_log_T2 %>% dplyr::select("anon_sample_ID", 1:13)



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


##### save both and load in 
write.csv(df_statsT2, "WHO_Serology/Limma_WHO_T2_PCVvPlacebo_anon.csv")


############### load in the limma files again and run the volcanos 
df<- read.csv("WHO_Serology/Limma_WHO_T2_PCVvPlacebo_anon.csv")




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

df$delabel<- gsub("IgG_", "" ,df$delabel)
df$delabel<- gsub("type.", "" ,df$delabel)


VP <- ggplot(data = df, aes(x = logFC, y = -log10(adj.P.Val), label = delabel)) + 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 9) +
  theme_classic() +
  ggtitle("T2 Placebo V PCV13") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", col = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  scale_fill_manual(values = mycolors) +
  geom_label_repel(size = 7, force = 50, max.overlaps = 25) +            # slightly bigger labels
  xlim(c(-2, 2)) +
  ylim(c(0, 2.5))

VP

ggsave(plot = VP, "Figures/Serology/Paper_final/VP_WHO_T2_PCVvPlac_PCV.png", dpi = 1200, height = 9, width = 9, units = "in")



######################## serology boxplots ################################## 


WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")
WHO_T0_T2
WHO_T0_T2<- WHO_T0_T2[,-c(1,2)]

WHO_T0_T2$Status<- factor(WHO_T0_T2$Status, levels = c( "PCV13 T0", "PCV13 T2", "Placebo T0", "Placebo T2")) 


#my_comps<- list(c("Placebo T0", "Placebo T2"), c("PCV13 T0", "PCV13 T2"), c("Placebo T0", "PCV13 T0"), c("Placebo T2", "PCV13 T0") )
my_comps_2<- list( c("PCV13 T0", "PCV13 T2"))
my_comps_3<- list( c("PCV13 T0", "PCV13 T2"), c("Placebo T0", "Placebo T2"))

WHO_T0_T2

colnames(WHO_T0_T2)

table(WHO_T0_T2$Status)



# Split data by species
iris_split <- split(WHO_T0_T2$IgG_type.6B, WHO_T0_T2$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 

WHO_T0_T2$IgG_type.14
################### boxplot Serotpye with cohort_time ######

serotype <- ggplot(WHO_T0_T2, aes(x = Status, y = IgG_type.19F)) +
  geom_boxplot(aes(fill = Status), width = 0.35, outlier.shape = NA,
               shape = 21, alpha = 0.4) +
  geom_jitter(aes(fill = Status), width = 0.25, shape = 21, alpha = 1, size = 4) +
  scale_fill_manual(values = c(
    "Placebo T0" = "#bbcfc7",
    "Placebo T2" = "#577a6c",
    "PCV13 T0"   = "#f3c6d9",
    "PCV13 T2"   = "#d32c75"
  )) +
  theme_classic() +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  ylab("") + #IgG ug/ml
  xlab("") +
  ggtitle("19F") +
  geom_vline(xintercept = 2.5, linetype = "dotted", color = "black") +
  geom_signif(
    comparisons = my_comps_2,
    map_signif_level = TRUE,
    test = "wilcox.test",
    textsize = 9
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  )


serotype

#c("Placebo T0", "Placebo T2"),
ggsave(plot = serotype, "Figures/Serology/Paper_final/Boxplots_paper/Serotype19F_box_paper.png", dpi = 1200,  height = 9,  width = 9)
######## extract just the legend 

serotype_leg <- serotype +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,          # force single row
      byrow = TRUE
    )
  )


library(cowplot)

legend_only <- get_legend(serotype_leg)
plot(legend_only)

ggsave(
  "Figures/Serology/Paper_final/Boxplots_paper/serotype_legend.png",
  legend_only,
  width = 10,
  height = 2,
  dpi = 1200
)

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
       y = "Serotype") +
  theme_bw() +
  scale_y_discrete(limits = rev(levels(factor(long_format$IgG_type))))+
  scale_size(range = c(0.1, 25)) +
  theme(
    plot.title   = element_text(size = 18, hjust = 0.5, face = "bold"), # title
    axis.title.x = element_text(size = 18, face = "bold"),              # x-axis title
    axis.title.y = element_text(size = 18, face = "bold"),              # y-axis title
    axis.text.x  = element_text(size = 16, face = "bold"),              # x tick labels
    axis.text.y  = element_text(size = 16, face = "bold"),              # y tick labels
    legend.title = element_text(size = 16, face = "bold"),              # legend title
    legend.text  = element_text(size = 16)               # legend labels
  )


map<- heat+ scale_color_viridis_c(option = "inferno")
map

ggsave(plot = map, "Figures/Serology/Paper_final/Bubble_ratio_PCV_placebo_new.png", dpi = 1200, height = 12, width =12, units = "in")




########## correlation of FC changes ##########################
library(corrplot)
library(ggcorrplot)

df_FC <- read.csv("WHO_Serology/Fold_change_df_anon_2.csv")
df_FC

colnames(df_FC)

FC_matrix<- df_FC %>% dplyr::select(patient_num , IgG_type, FC_T2)


FC_matrix
FC_matrix_PCV<- FC_matrix 
#FC_matrix_PCV<- FC_matrix[FC_matrix$Treatment== "PCV13", ]

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

colnames(FC_matrix_pivot)

FC_matrix_pivot_2<- FC_matrix_pivot %>% column_to_rownames(var="patient_num")

#FC_matrix_pivot_2<- FC_matrix_pivot_2 %>% select(!Treatment  )

FC_matrix_pivot_2
FC_matrix_pivot_2<- as.matrix(FC_matrix_pivot_2)


FC_matrix_pivot_2


# Create the correlation matrix
cor_matrix <- cor(FC_matrix_pivot_2, use = "pairwise.complete.obs")
cor_matrix


# Plot the correlation matrix using ggcorrplot


corrplot(cor_matrix, type = 'lower', order = 'hclust', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45, col = rev(COL2('RdBu', 10)))

corrplot(
  cor_matrix,
  type = 'lower',
  order = 'hclust',
  tl.col = 'black',       # text label color
  tl.srt = 45,            # text rotation
  tl.cex = 1.5,           # text label size
  cl.cex = 1.2,           # color legend size
  number.cex = 1.2,       # correlation number size (if shown)
  cl.ratio = 0.2,         
  col = rev(COL2('RdBu', 10))
)

# Create a PDF file to save the plot
pdf("Figures/Serology/Paper_final/FC_corrplot_update.pdf", width = 10, height = 10)  # Specify the file name and dimensions in inches

corrplot(
  cor_matrix,
  type = 'lower',
  order = 'hclust',
  tl.col = 'black',       # text label color
  tl.srt = 45,            # text rotation
  tl.cex = 1.5,           # text label size
  cl.cex = 1.2,           # color legend size
  number.cex = 1.2,       # correlation number size (if shown)
  cl.ratio = 0.2,         
  col = rev(COL2('RdBu', 10))
)
# Close the PDF device to save the plot
dev.off()


