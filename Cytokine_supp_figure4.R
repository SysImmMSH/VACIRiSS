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
my_comps <- list( c("PCV13 T0", "PCV13 T2"), c("PCV13 T1", "PCV13 T2"), c("Placebo T0", "Placebo T1"), c("Placebo T0", "Placebo T2"), c("PCV13 T0", "Healthy"),  c("PCV13 T1", "Healthy"), c("Placebo T0", "Healthy"), c("Placebo T1", "Healthy"))


my_comps <- list( c("PCV13 T0", "PCV13 T1"), c("PCV13 T0", "PCV13 T2"))


my_comps <- list( c("PCV13 T0", "Healthy"),  c("PCV13 T1", "Healthy"), c("PCV13 T2", "Healthy"),c("Placebo T0", "Placebo T2"), c("Placebo T0", "Healthy"), c("Placebo T1", "Healthy"), c("Placebo T2", "Healthy") )

my_comps_healthy <- list( c("PCV13 T0", "Healthy"),   c("PCV13 T2", "Healthy"), c("Placebo T0", "Healthy"),  c("Placebo T1", "Healthy"),  c("Placebo T2", "Healthy"))

colnames(cyto_health)



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



cyto_plot<-ggplot(cyto_health, aes(x= Status, y= IFNy)) +
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
  ylab("pg/ml")+
  xlab("")+
  ggtitle("IFNy")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black") +
  geom_vline(xintercept = 1.5, linetype= "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_healthy, method= "wilcox.test",  label = "p.signif", hide.ns = TRUE) +
  scale_y_continuous(trans = "log1p")

#geom_signif(comparisons = my_comps2, map_signif_level = TRUE, test = "wilcox.test")

cyto_plot

#c("Placebo T0", "Placebo T2"),
ggsave(plot = cyto_plot, "Figures/Cytokine/Healthy_Timepoint/CXCL10_box_healthy.png", dpi = 300,  height = 7,  width = 7)

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


###################################  responder v non 


# Remove rows with NA in responder
colnames(cyto_health)
cyto_resp <- cyto_health %>% filter(!is.na(Responder_final))
cyto_resp


cyto_resp$Responder_Timepoint<- paste(cyto_resp$Timepoint, cyto_resp$Responder_final, sep = "_")

cyto_resp$Responder_Timepoint
cyto_resp$Status

colnames(cyto_resp)

# Split data by species
iris_split <- split(cyto_resp$sCD40L, cyto_resp$Responder_Timepoint)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 


colnames(cyto_resp)

colnames(cyto_resp)

table(cyto_resp$Responder_Timepoint)

t_test_results <- compare_means(sCD40L ~ Responder_Timepoint, data = cyto_resp, method = "wilcox.test")
t_test_results


significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results

Comp_resp_all<- list(c("T0_Yes", "T0_No"), c("T1_Yes", "T1_No"), c("T2_Yes", "T2_No"), c("T0_Yes", "T1_Yes"), c("T0_Yes", "T2_Yes"), c("T1_Yes", "T2_Yes"),  c("T0_No", "T1_No"), c("T0_No", "T2_No"), c("T1_No", "T2_No"))

Comp_resp<- list(c("T2_Yes", "T2_No"), c("T0_Yes", "T2_Yes"), c("T1_Yes", "T2_Yes"), c("T0_No", "T2_No"))

Comp_resp<- list(c("T2_Yes", "T2_No"), c("T0_Yes", "T2_Yes"), c("T1_Yes", "T2_Yes"), c("T0_No", "T2_No"))

Comp_resp<- list( c("T0_No", "T2_No"), c("T1_No", "T2_No"), c("T2_Yes", "T2_No"))


Comp_resp<- list( c("T2_Yes", "T2_No"),  c("T0_Yes", "T2_Yes"))
Comp_resp<- list( c("T0_Yes", "T2_Yes"))
Comp_resp<- list( c("T0_No", "T2_No"))

Comp_resp<- list(  c("T2_No", "T2_Yes"), c("T1_Yes", "T2_Yes"), c("T0_Yes", "T2_Yes"))



cyto_resp

flow_plot<-ggplot(cyto_resp, aes(x= Responder_Timepoint, y= `sCD40L`)) +
  geom_boxplot(aes(fill= Responder_Timepoint), width= 0.35, outlier.shape = NA, alpha = 0.4)+
  #geom_hline(yintercept = 54:68, linetype= "dashed", color = "grey",  alpha = 0.6)+
  #geom_jitter(aes(fill= Status, shape= Responder), width=0.25, alpha=1, size = 3) +
  geom_jitter(aes(fill= Responder_Timepoint),  shape= 21, width=0.25, alpha=1, size = 3) +
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
  #scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  #annotation_ticks(sides = "l", type = "minor", outside = TRUE)+
  ylab("pg/ml")+
  xlab("")+
  ggtitle("sCD40L")+
  #ylim(0, 100) +
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black")  +
  #geom_vline(xintercept = 1.5, linetype= "dotted", color = "black")  +
  #geom_hline(yintercept = 54:68, linetype= "dashed", color = "red")  +
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = Comp_resp, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE) + #label.y = c(95, 90, 85, 80, 75, 40) 
  scale_y_continuous(trans = "log1p")

flow_plot

ggsave(plot = flow_plot, "Figures/Responder/Cytokine/sCD40L.png", dpi = 300,  height = 6,  width = 6)
