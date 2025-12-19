############## cytokine. Vaciriss biopaper ###
## comparing healthy controls from immerse to vacciriss ####
#### volcano plot #####
## boxplot ###


library(tidyverse)
library(ggplot2)
library(readr)
library(limma)
library(ggrepel)
library(ggpubr)
library(readxl)
library(rstatix)
library(pheatmap)
library(ComplexHeatmap)
library(tidyverse)
library(ggbeeswarm)
library(wesanderson)
library(plyr)
library(dplyr)
library(ggpubr)
library(cowplot)
library(rlang)
#install.packages("ggtext")
library(ggtext)
library("factoextra")
library(umap)

#### load in the data 

getwd()

cyto_health<- read.csv("paper/final anonymised dataframes/cytokine/cytokine_anon.csv")


cyto_health
cyto_health$Status 


###### plot the healthy 




my_comps_healthy <- list( c("PCV13 T0", "Healthy"),  c("PCV13 T1", "Healthy"),  c("PCV13 T2", "Healthy"), c("Placebo T0", "Healthy"),  c("Placebo T1", "Healthy"),  c("Placebo T2", "Healthy"))
my_comps <- list( c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"), c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("PCV13 T0", "Healthy"),  
                  c("PCV13 T1", "Healthy"), c("Placebo T0", "Healthy"), c("Placebo T1", "Healthy"))


my_comps <- list( c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"))


my_comps <- list( c("PCV13 T0", "Healthy"),  c("PCV13 T1", "Healthy"), c("PCV13 T2", "Healthy"),c("Placebo T0", "Placebo T2"), c("Placebo T0", "Healthy"), c("Placebo T1", "Healthy"), c("Placebo T2", "Healthy") )

my_comps_healthy <- list( c("PCV13 T0", "Healthy"),   c("PCV13 T2", "Healthy"), c("Placebo T0", "Healthy"),  c("Placebo T1", "Healthy"),  c("Placebo T2", "Healthy"))

colnames(cyto_health)

my_comps_IL10 <- list( c("PCV13 T0", "PCV13 T2"), c("PCV13 T0", "Healthy"),  c("PCV13 T1", "Healthy"), c("Placebo T0", "Placebo T2"), c("Placebo T0", "Healthy"), c("Placebo T1", "Healthy") )




# Split data by species
iris_split <- split(cyto_health$IFNy, cyto_health$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 



colnames(cyto_health)

t_test_results <- compare_means(IFNy ~ Status, data = cyto_health, method = "wilcox.test")
t_test_results


significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results


colnames(cyto_health)

cyto_plot<-ggplot(cyto_health, aes(x= Status, y= sCD40L)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Status), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("Placebo T0" = "#bbcfc7",
                              "Placebo T1" = "#85a89a",
                              "Placebo T2" ="#577a6c",
                              "Placebo T3" ="#30443c",
                              "PCV13 T0" = "#f3c6d9", 
                              "PCV13 T1" = "#e379a7",
                              "PCV13 T2" = "#d32c75",
                              "PCV13 T3" = "#ac2460",
                              "Healthy" = "grey"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("")+ #pg/ml
  xlab("")+
  ggtitle("sCD40L")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  ) +
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black") +
  geom_vline(xintercept = 1.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_healthy, method= "wilcox.test",  label = "p.signif", hide.ns = TRUE, size= 9) +
  scale_y_continuous(trans = "log1p")

#geom_signif(comparisons = my_comps2, map_signif_level = TRUE, test = "wilcox.test")

cyto_plot

#c("Placebo T0", "Placebo T2"),
ggsave(plot = cyto_plot, "Figures/Cytokine/Healthy_Timepoint/Paper/sCD40L_box_healthy.png", dpi = 1200,  height = 9,  width = 9)

###### stats 

colnames(cyto_health)

summary <- cyto_health %>%
  group_by(Status) %>%
  dplyr::summarise(
    median_percentage = median(IFNy, na.rm = TRUE),
    Q1 = quantile(IFNy, 0.25, na.rm = TRUE),
    Q3 = quantile(IFNy , 0.75, na.rm = TRUE),
    sample_count=n()
  )

summary
