############ Vaciriss ##
##### boxplots for flow bring in the responder status as well ###
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(readr)
#install.packages("ggpubr")
library(ggpubr)
library(dplyr)
library(stringr)
library(rstatix)
#install.packages("ggprism")
library(ggprism)

getwd()

Bcell_df_HTC<- read.csv("paper/final anonymised dataframes/flow/Bcell_flow_anon.csv")

## filter it 


Bcell_df_HTC_filter<- Bcell_df_HTC %>% filter(CD19.Count >1000)

Bcell_df_HTC_filter

colnames(Bcell_df_HTC_filter)

# Split data by species
iris_split <- split(Bcell_df_HTC_filter$B1, Bcell_df_HTC_filter$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if greater than 0.05 then it is not normal 



############## to do them individually 
t_test_results <- compare_means(B1 ~ Status, data = Bcell_df_HTC_filter, method = "wilcox.test")
t_test_results

significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results

#### check the stats are right 
my_comps_5<- list(c("Healthy", "PCV13 T0"), c("Healthy", "PCV13 T1"), c("Healthy", "PCV13 T2"),  c("Healthy", "Placebo T0"), c("Healthy", "Placebo T1"), c("Healthy", "Placebo T2")) ## healthy v everything

Bcell_df_HTC_filter
colnames(Bcell_df_HTC_filter)

flow_plot<-ggplot(Bcell_df_HTC_filter, aes(x= Status, y= `Naive.mature`)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, alpha = 0.4)+
  geom_jitter(aes(fill= Status),  shape= 21, width=0.25, alpha=1, size = 3) +
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
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black")
  ) +
  ylab("% of CD19+")+
  xlab("")+
  ggtitle("IgG Memory")+
  #ylim(0, 50) +
  theme(plot.title = element_text(hjust=0.5))+
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_5, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE, label.y = c(30, 40))  #+ #label.y = c(95, 90, 85, 80, 75, 40) 



flow_plot

ggsave(plot = flow_plot, "Figures/Flow/updated_flow/Healthy_boxplots/Bcells/Paper/IgG_memory_total.png", dpi = 300,  height = 6,  width = 6)

######## summarise the columns 
# Summarize the mean of a specific column (e.g., CD19.Count) grouped by Status
summary <- Bcell_df_HTC_filter %>%
  group_by(Status) %>%
  dplyr::summarise(
    median_percentage = median(Naive.mature, na.rm = TRUE),
    Q1 = quantile(Naive.mature, 0.25, na.rm = TRUE),
    Q3 = quantile(Naive.mature, 0.75, na.rm = TRUE),
    sample_count=n()
  )

summary


################################################################
########################## T cells ##############################
################################################################

Tcell_df_HTC<- read.csv("paper/final anonymised dataframes/flow/Tcell_flow_anon.csv")
nrow(Tcell_df_HTC)

Tcell_df_HTC_filter<- Tcell_df_HTC %>% filter(CD3.Count >1000)
nrow(Tcell_df_HTC_filter)

colnames(Tcell_df_HTC_filter)

# Split data by species
iris_split <- split(Tcell_df_HTC_filter$Th1, Tcell_df_HTC_filter$Status)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if under 0.05 then it is not normal 


Tcell_df_HTC_filter
t_test_results <- compare_means(Th1 ~ Status, data = Tcell_df_HTC_filter, method = "t.test")
t_test_results

significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results



#### check the stats are right 
colnames(Tcell_df_HTC_filter)
my_comps_5<- list(c("Healthy", "PCV13 T0"), c("Healthy", "PCV13 T1"), c("Healthy", "PCV13 T2"),  c("Healthy", "Placebo T0"), c("Healthy", "Placebo T1"), c("Healthy", "Placebo T2")) ## healthy v everything

my_comps_5<- list(c("Healthy", "PCV13 T0"), c("Healthy", "PCV13 T1"), c("Healthy", "PCV13 T2"),  c("Healthy", "Placebo T0"), c("Healthy", "Placebo T1"), c("Healthy", "Placebo T2")) ## healthy v everything


flow_plot<-ggplot(Tcell_df_HTC_filter, aes(x= Status, y= `Th1`)) +
  geom_boxplot(aes(fill= Status), width= 0.35, outlier.shape = NA, alpha = 0.4)+
  geom_jitter(aes(fill= Status),  shape= 21, width=0.25, alpha=1, size = 3) +
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
  theme(
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black")
  ) +
  ylab("% of CD8")+
  xlab("")+
  ggtitle("CD8 Effector Memory")+
  theme(plot.title = element_text(hjust=0.5))+
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_5, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE )# + #label.y = c(95, 90, 85, 80, 75, 40) 

flow_plot

ggsave(plot = flow_plot, "Figures/Flow/updated_flow/Healthy_boxplots/Tcells/CD4 counts/CD8.effector.memory.png", dpi = 300,  height = 6,  width = 6)


colnames(Tcell_df_HTC_filter)
######## summarise the columns 
# Summarize the mean of a specific column (e.g., CD19.Count) grouped by Status
summary <- Tcell_df_HTC_filter %>%
  group_by(Status) %>%
  dplyr::summarise(
    median_percentage = median(Th1, na.rm = TRUE),
    Q1 = quantile(Th1, 0.25, na.rm = TRUE),
    Q3 = quantile(Th1, 0.75, na.rm = TRUE),
    sample_count=n()
  )

summary



###################### responder analysis ##################


Bcell_df_HTC_filter 

Bcell_df_HTC_filter$Responder_timepoint<- paste(Bcell_df_HTC_filter$Timepoint, Bcell_df_HTC_filter$Responder_final, sep = "_")

colnames(Bcell_df_HTC_filter)
Bcell_df_HTC_filter
# Remove rows with NA in responder
df_clean <- Bcell_df_HTC_filter %>% filter(!is.na(Responder_final))
df_clean


colnames(df_clean)
# Split data by species
iris_split <- split(df_clean$IgG_memory_total, df_clean$Responder_timepoint)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 



colnames(df_clean)
t_test_results <- compare_means(IgG_memory_total ~ Responder_timepoint, data = df_clean, method = "t.test")
t_test_results


significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results

Comp_resp_1<- list( c("T1_No", "T2_No"))
Comp_resp_1<- list( c("T0_Yes", "T0_No"))



flow_plot<-ggplot(df_clean, aes(x= Responder_timepoint, y= `IgG_memory_total`)) +
  geom_boxplot(aes(fill= Responder_timepoint), width= 0.35, outlier.shape = NA, alpha = 0.6)+
  #geom_hline(yintercept = 54:68, linetype= "dashed", color = "grey",  alpha = 0.6)+
  #geom_jitter(aes(fill= Status, shape= Responder), width=0.25, alpha=1, size = 3) +
  geom_jitter(aes(fill= Responder_timepoint),  shape= 21, width=0.25, alpha=1, size = 3) +
  #scale_shape_manual(values = c("Yes" = 21, "No" = 24))+
  scale_fill_manual(values= c(
    "T0_No" = "#02dbc7", 
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
  #annotation_logticks(sides = "l")+
  #annotation_ticks(sides = "l", type = "minor", outside = TRUE)+
  ylab("% of CD19")+
  xlab("")+
  ggtitle("IgG Memory")+
  #ylim(0, 20) +
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black")  +
#geom_vline(xintercept = 1.5, linetype= "dotted", color = "black")  +
#geom_hline(yintercept = 54:68, linetype= "dashed", color = "red")  +
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = Comp_resp_1, method= "t.test",  label = "p.signif", hide.ns = FALSE)# +
#stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = Comp_resp_2, method= "t.test",  label = "p.signif", hide.ns = FALSE) + #+#label.y = c(95, 90, 85, 80, 75, 40) 
#scale_y_continuous(trans = "log1p", limits = c(0,10)) #limits = c(0,100)

flow_plot

ggsave(plot = flow_plot, "Figures/Responder/Flow/IgG.Memory_responder.png", dpi = 300,  height = 6,  width = 6)
######### summary ####

# Summarize the mean of a specific column (e.g., CD19.Count) grouped by Status
summary <- df_clean %>%
  group_by(Responder_timepoint) %>%
  dplyr::summarise(
    median_percentage = median(Unswitched.memory, na.rm = TRUE),
    Q1 = quantile(Unswitched.memory, 0.25, na.rm = TRUE),
    Q3 = quantile(Unswitched.memory, 0.75, na.rm = TRUE),
    sample_count=n()
  )

summary




########. Repeat for T cells ######

Tcell_df_HTC_filter



Tcell_df_HTC_filter$Responder_timepoint<- paste(Tcell_df_HTC_filter$Timepoint, Tcell_df_HTC_filter$Responder_final, sep = "_")


# Remove rows with NA in a specific column
df_clean_T <- Tcell_df_HTC_filter %>% filter(!is.na(Responder_final))
df_clean_T




table(df_clean_T$Responder_timepoint)

colnames(df_clean_T)
# Split data by species
iris_split <- split(df_clean_T$CD8.Naive, df_clean_T$Responder_timepoint)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 



colnames(df_clean_T)
t_test_results <- compare_means(CD8.Naive ~ Responder_timepoint, data = df_clean_T, method = "t.test")
t_test_results


significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results

Comp_resp_1<- list( c("T0_No", "T0_Yes"))
Comp_resp_3<- list( c("T0_No", "T2_No"))

Comp_resp_2<- list(  c("T0_No", "T1_No"), c("T0_No", "T2_No"))

Comp_resp_1<- list( c("T0_Yes", "T2_Yes"), c("T2_Yes", "T2_No"))
Comp_resp_3<- list( c("T0_No", "T1_No"))
Comp_resp_3<- list( c("T0_Yes", "T0_No"), c("T1_Yes", "T1_No"), c("T2_Yes", "T2_No"))


Comp_resp_1<- list( c("T0_Yes", "T2_Yes"), c("T2_Yes", "T2_No"))

Comp_resp_1<- list( c("T2_Yes", "T0_Yes"))
Comp_resp_2<- list( c("T0_Yes", "T1_Yes"))

Comp_resp_1<- list( c("T2_No", "T0_No"))
Comp_resp_1<- list( c("T2_No", "T0_No"))



flow_plot<-ggplot(df_clean_T, aes(x= Responder_timepoint, y= `CD8.Naive`)) +
  geom_boxplot(aes(fill= Responder_timepoint), width= 0.35, outlier.shape = NA, alpha = 0.6)+
  #geom_hline(yintercept = 54:68, linetype= "dashed", color = "grey",  alpha = 0.6)+
  #geom_jitter(aes(fill= Status, shape= Responder), width=0.25, alpha=1, size = 3) +
  geom_jitter(aes(fill= Responder_timepoint),  shape= 21, width=0.25, alpha=1, size = 3) +
  #scale_shape_manual(values = c("Yes" = 21, "No" = 24))+
  scale_fill_manual(values= c(
    "T0_No" = "#02dbc7", 
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
  #annotation_logticks(sides = "l")+
  #annotation_ticks(sides = "l", type = "minor", outside = TRUE)+
  ylab("% of CD8")+
  xlab("")+
  ggtitle("CD8 Naive")+
  #ylim(0, 20) +
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
#geom_vline(xintercept = 1.5, linetype= "dotted", color = "black")  +
#geom_hline(yintercept = 54:68, linetype= "dashed", color = "red")  +
#stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = Comp_resp_1, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE) #+
stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = Comp_resp_2, method= "t.test",  label = "p.signif", hide.ns = FALSE)# + #+#label.y = c(95, 90, 85, 80, 75, 40) 
#scale_y_continuous(trans = "log1p", limits = c(60,99)) #limits = c(0,100)

flow_plot

ggsave(plot = flow_plot, "Figures/Responder/Flow/CD8.Naive_responder.png", dpi = 300,  height = 6,  width = 6)

#######

# Summarize the mean of a specific column (e.g., CD19.Count) grouped by Status
summary <- df_clean_T %>%
  group_by(Responder_timepoint) %>%
  dplyr::summarise(
    median_percentage = median(CD8.Naive, na.rm = TRUE),
    Q1 = quantile(CD8.Naive, 0.25, na.rm = TRUE),
    Q3 = quantile(CD8.Naive, 0.75, na.rm = TRUE),
    sample_count=n()
  )

summary


table(df_clean_T$Responder_timepoint)
