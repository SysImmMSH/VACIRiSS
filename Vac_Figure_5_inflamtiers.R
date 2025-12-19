###### 
library(stringr)
library(dplyr)
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


########

df<- read.csv("Transcript/new_analysis/vsd_batch_corrected.csv")

df

metadata <- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata.csv") #400
metadata

metadata <- metadata %>% dplyr::select(anon_sample_ID, Timepoint, Treatment, Cohort_time, patient_num)
metadata


inflm_meta<- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata_withinflamtiers.csv")
inflm_meta

####### Ig and cytotoxic gene levels with inflam score 

# Vector of gene names you want to keep
selected_genes <- c("IGHM", "IGHA2", "IGHG2", "IGHA1", "IGHG1", "IGHG3","NCAM1", "PRF1", "GNLY")

# Filter rows where gene_name is in selected_genes
Ig_genes <- df[df$gene_name %in% selected_genes, ]

colnames(Ig_genes)
Ig_genes
Ig_genes<- Ig_genes[,-1]

#### flip this to merge it with the meta
# Step 1: Set gene names as rownames
rownames(Ig_genes) <- Ig_genes$gene_name

# Step 2: Drop the 'gene_name' column to keep only numeric values
Ig_genes <- Ig_genes[, -which(names(Ig_genes) == "gene_name")]

# Step 3: Transpose the data
df_t <- as.data.frame(t(Ig_genes))

df_t
# Step 4: Optional - make gene names proper column names
#colnames(df_t) <- rownames(df)

# Step 5: Add sample names as a column if you want
df_t$anon_sample_ID <- rownames(df_t)

# Step 6: Reorder to put sample column first
#df_t <- df_t[, c(ncol(df_t), 1:(ncol(df_t) - 1))]

# View result
head(df_t)


metadata_filt<- inflm_meta %>% dplyr::select(anon_sample_ID, Timepoint,Treatment, Cohort_time,  Responder_final, patient_num, inflam.gs)

IgG_meta<- merge(df_t, metadata_filt, by= "anon_sample_ID")

IgG_meta

IgG_meta
colnames(IgG_meta)

##### keep only paired T0 and T1 samples


IgG_meta_T0_T1<- IgG_meta[IgG_meta$Timepoint == "T0"| IgG_meta$Timepoint == "T1",]

table(IgG_meta_T0_T1$Timepoint)

IgG_meta_T0_T1

# Count how many times each patient appears
paired_patients <- IgG_meta_T0_T1 |>
  dplyr::group_by(patient_num) |>
  dplyr::summarise(n_timepoints = n_distinct(Timepoint)) |>
  dplyr::filter(n_timepoints >= 2)

paired_patients
# Filter original dataframe to keep only those patients
df_paired <- IgG_meta_T0_T1[IgG_meta_T0_T1$patient_num %in% paired_patients$patient_num, ]

df_paired
###### log 2 

df_paired_log<- df_paired
df_paired

df_paired_log[c("IGHA2", "IGHG2", "IGHA1", "IGHG1", "IGHG3", "IGHM", "NCAM1", "PRF1", "GNLY")] <- log2(df_paired_log[c("IGHA2", "IGHG2", "IGHA1", "IGHG1", "IGHG3", "IGHM", "NCAM1", "PRF1", "GNLY")] + 1)


my_comps <- list( c("PCV13_T0", "PCV13_T1"), c("Placebo_T0", "Placebo_T1"))
######### significance and normality 

# Split data by species
iris_split <- split(df_paired_log$NCAM1, df_paired_log$Cohort_time)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 



colnames(cyto_health)

t_test_results <- compare_means(NCAM1 ~ Cohort_time, data = df_paired_log, method = "t.test")
t_test_results


significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results

df_paired_log

#dev.off()
colnames(df_paired_log)
####### ggplot 

cyto_plot <- ggplot(df_paired_log, aes(x = Cohort_time, y = IGHA2)) +
  geom_boxplot(aes(fill = Cohort_time), width = 0.35, outlier.shape = NA, shape = 21, alpha = 0.4) +
  geom_jitter(aes(fill = Cohort_time), width = 0.25, shape = 21, alpha = 1, size = 4) +
  #geom_point(aes(fill = Cohort_time), width = 0.25, shape = 21, alpha = 1, size = 4) +
  geom_line(aes(group = patient_num), color = "gray60", linewidth = 0.5, alpha = 0.3) +  # ← Add this line for connecting dots
  scale_fill_manual(values = c(
    "Placebo_T0" = "#bbcfc7",
    "Placebo_T1" = "#85a89a",
    "Placebo_T2" = "#577a6c",
    "Placebo_T3" = "#30443c",
    "PCV13_T0" = "#f3c6d9",
    "PCV13_T1" = "#e379a7",
    "PCV13_T2" = "#d32c75",
    "PCV13_T3" = "#ac2460",
    "Healthy" = "grey"
  )) +
  theme_classic() +
  annotation_logticks(sides = "l") +
  ylab("log2 normalised counts") +
  xlab("") +
  ggtitle("IGHA2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 2.5, linetype = "dotted", color = "black") +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps, method= "t.test",  label = "p.signif", hide.ns = TRUE) #+ p.stat for the nu,bers

#scale_y_continuous(trans = "log1p")

cyto_plot

#c("Placebo T0", "Placebo T2"),
#ggsave(plot = cyto_plot, "Figures/Ravichan_tests/Gene_expression_boxplots/NCAM1_boxplot.png", dpi = 300,  height = 7,  width = 7)


###### add the inflamation tiers 


df_paired_log$Time_inf<- paste(df_paired_log$Timepoint, df_paired_log$inflam.gs, sep = "_")

####### this one for 
df_paired_log$Time_inf<- factor(df_paired_log$Time_inf, levels = c( "T0_low", "T1_low", "T0_high", "T1_high")) 



df_paired_log_PCV<- df_paired_log[df_paired_log$Treatment == "PCV13",]

my_comparisons <- list(
  c("T0_low", "T1_low"),
  c("T0_high", "T1_high"),
  c("T0_low", "T0_high"),
  c("T1_low", "T1_high"))


my_comparisons_2 <- list(
  c("T0_low", "T1_low"))



df_paired_log
colnames(df_paired_log_PCV)


# Split data by species
iris_split <- split(df_paired_log_PCV$PRF1, df_paired_log_PCV$Time_inf)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if less than 0.05 then it is not normal 

colnames(df_paired_log_PCV)

flow_plot <- ggplot(df_paired_log_PCV, aes(x = Time_inf, y = PRF1)) +
  geom_boxplot(aes(fill = Time_inf), width = 0.35,
               outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill = Time_inf), shape = 21,
              width = 0.25, alpha = 1, size = 3) +
  scale_fill_manual(
    name = "Inflammation tier",
    values = c(
      "T0_low"  = "#FFEB3B",
      "T1_low"  = "#FFEB3B",
      "T0_high" = "red",
      "T1_high" = "red"
    ),
    labels = c(
      "T0_low"  = "T0 (Low)",
      "T1_low"  = "T1 (Low)",
      "T0_high" = "T0 (High)",
      "T1_high" = "T1 (High)"
    )
  ) +
  scale_x_discrete(labels = c(
    "T0_low"  = "T0 (Low)",
    "T1_low"  = "T1 (Low)",
    "T0_high" = "T0 (High)",
    "T1_high" = "T1 (High)"
  )) +
  theme_classic() +
  ylab("") + #log2 normalised counts
  xlab("") +
  ggtitle("PRF1") +
  
  geom_vline(xintercept = 2.5, linetype = "dotted", color = "black") +
  
  stat_compare_means(
    aes(label = ..p.signif..),
    comparisons = my_comparisons,
    method = "wilcox.test",
    hide.ns = FALSE,
    size = 9
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
#+ 
#stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comparisons_2, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE )# + #label.y = c(95, 90, 85, 80, 75, 40) 


flow_plot

ggsave(plot = flow_plot, "Figures/Modules/paper/Infltiers_timepoint_PRF1.png", dpi = 1200,  height = 9,  width = 9)


#### get the legend by itself 


flow_plot_leg <- flow_plot +
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

legend_only <- get_legend(flow_plot_leg)
plot(legend_only)

ggsave(
  "Figures/Modules/paper/inflam_legend.png",
  legend_only,
  width = 10,
  height = 2,
  dpi = 1200
)


##################################################################
######################## Heatmap of IgG ################################


head(df_paired_log)

colnames(df_paired_log)



rownames(df_paired_log)<- NULL
Ig_matrix<- df_paired_log %>% dplyr::select(anon_sample_ID, IGHA2, IGHG2, IGHA1, IGHG1, IGHG3, IGHM) %>% 
  column_to_rownames("anon_sample_ID")


Ig_matrix<- as.matrix(Ig_matrix)
Ig_matrix

Ig_matrix_fixed<- t(Ig_matrix)

Ig_meta<- df_paired_log %>% dplyr::select(anon_sample_ID, Treatment, inflam.gs) %>% 
  column_to_rownames("anon_sample_ID")


table(Ig_meta$Treatment)
Ig_meta



library(ComplexHeatmap)


dim(Ig_meta)
dim(Ig_matrix)
#coldata<- Ig_meta
#coldata
Ig_meta_treatment<- Ig_meta %>% dplyr::select("Treatment")

Ig_meta_treatment
#colnames(coldata)[which(names(coldata) == "Cohort_time")] <- "Timepoint by Treatment"
# Column annotation
col_ha <- HeatmapAnnotation(df = Ig_meta_treatment)
col_ha

library(circlize)

# Define your custom color palette
custom_colors <- colorRampPalette(c("#163f76", "#4575B4", "#FFFFF2", "#ec763f", "#bf2a1b"))(100)


# Construct metadata
coldata <- data.frame(Ig_meta_treatment)
rownames(Ig_meta_treatment)

ha_colors <- list(
  Treatment = c("Placebo" = "#30443c", "PCV13" = "#ac2460")
)


#p this palette to the range of expression values
range_vals <- range(Ig_matrix_fixed, na.rm = TRUE)
print(range_vals)

col_fun <- colorRamp2(
  seq(from = range_vals[1], to = range_vals[2], length.out = 100),
  custom_colors
)


col_ha <- HeatmapAnnotation(
  df = coldata,
  col = ha_colors
)


dim(Ig_matrix_fixed)

# Split the heatmap by Timepoint
H <- Heatmap(
  Ig_matrix_fixed,
  top_annotation = col_ha,
  column_split = Ig_meta_treatment$Treatment,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  name = "Expression",
  col = col_fun,
  show_heatmap_legend = TRUE, 
  column_title = "T0 Immunoglobulin gene expression",
  row_names_side = "right",
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  row_names_gp = gpar(fontsize = 15),
  width = unit(5, "in"),
  height = unit(2, "in")
)

draw(H)

pdf("Figures/Modules/paper/heatmap_IgGenes_treatmentT0.pdf", width = 8, height = 10)  # set desired size
draw(H)  # H is your Heatmap object
dev.off()


lgd=Legend(col_fun = col_fun, "Expression")

png("Figures/Modules/paper/heatmap_IgGenes_treatmentT0.png", width = 10, height = 3, units = "in" ,res =  1200)  # set desired size
draw(H)  # H is your Heatmap object
dev.off()

#### try and

###### now repeat for the inflam in PCV T0

rownames(df_paired_log)<- NULL
df_paired_log
df_paired_log_PCV_T0<- df_paired_log[df_paired_log$Cohort_time == "PCV13_T0",]
rownames(df_paired_log_PCV_T0)<- NULL

Ig_matrix_T0<- df_paired_log_PCV_T0 %>% dplyr::select(anon_sample_ID, IGHA2, IGHG2, IGHA1, IGHG1, IGHG3, IGHM) %>% 
  column_to_rownames("anon_sample_ID")

Ig_matrix_T0<- as.matrix(Ig_matrix_T0)
Ig_matrix_T0

Ig_matrix_T0_fixed <- t(Ig_matrix_T0)


Ig_meta_T0<- df_paired_log[df_paired_log$Cohort_time == "PCV13_T0",]

Ig_meta_T0

rownames(Ig_meta_T0)<- NULL

Ig_meta_T0<- Ig_meta_T0 %>% dplyr::select(anon_sample_ID, Treatment, inflam.gs) %>% 
  column_to_rownames("anon_sample_ID")

Ig_meta_inflam<- Ig_meta_T0 %>% dplyr::select("inflam.gs")

Ig_meta_inflam
#colnames(coldata)[which(names(coldata) == "Cohort_time")] <- "Timepoint by Treatment"
# Column annotation
col_ha <- HeatmapAnnotation(df = Ig_meta_inflam)
col_ha

library(circlize)


dim(Ig_matrix_T0_fixed)
Ig_matrix_T0_fixed
dim(Ig_meta_inflam)
Ig_matrix_T0_fixed
Ig_matrix_T0_fixed_t<- t(Ig_matrix_T0_fixed)
Ig_meta_inflam

dim(Ig_matrix_T0_fixed_t)
dim(Ig_matrix_T0_fixed)
dim(Ig_matrix_T0)  # should be 6 65
dim(Ig_matrix_T0_fixed) # should be 6 65


# Check:
ncol(Ig_matrix_T0_fixed) == nrow(Ig_meta_inflam)  # Should be TRUE

# And:
all(colnames(Ig_matrix_T0_fixed) == rownames(Ig_meta_inflam))  # Should be TRUE


# Define your custom color palette
custom_colors <- colorRampPalette(c("#163f76", "#4575B4", "#FFFFF2", "#ec763f", "#bf2a1b"))(100)


# Construct metadata
coldata <- data.frame(Ig_meta_inflam)
rownames(Ig_meta_inflam)

ha_colors <- list(inflam.gs = c("low" = "#FFEB3B", "high" = "red"))

col_ha <- HeatmapAnnotation(
  df = Ig_meta_inflam,
  col = ha_colors
)



#p this palette to the range of expression values
range_vals <- range(Ig_matrix_T0_fixed_t, na.rm = TRUE)
print(range_vals)

col_fun <- colorRamp2(
  seq(from = range_vals[1], to = range_vals[2], length.out = 100),
  custom_colors
)


# Make a factor for the split and set its name to change the legend title
inflam_factor <- factor(Ig_meta_inflam$inflam.gs)
attr(inflam_factor, "name") <- "Inflammatory tier"  # this sets the legend title

# Split the heatmap by inflam 

H <- Heatmap(
  Ig_matrix_T0_fixed,
  top_annotation = col_ha,
  column_split = inflam_factor,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  name = "Expression",
  col = col_fun,
  show_heatmap_legend = TRUE, 
  column_title = "Inflam Tier Gene Expression",
  row_names_side = "right",
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  row_names_gp = gpar(fontsize = 15),
  width = unit(5, "in"),
  height = unit(2, "in")
)

draw(H)

png("Figures/Modules/paper/heatmap_IgGenes_inflammtiers.png", width = 10, height = 3, units = "in" , res =  1200)  # set desired size
draw(H)  # H is your Heatmap object
dev.off()

####### cyotoxic gene module 


gene_list<- "TSPOAP1
NCAM1
GNLY
SLC25A10
SIGLEC17P
TRDC
KLRC1
SLFN13
KLRF1
ADAMTS1
NMUR1
CLIC3
IL2RB
CRIM1
CD38
LINC00299
TRGC1
PRR34
ARVCF
FEZ1
SPON2
MLC1
LPAL2
ERBB2
HDGFL3
CD247
PRF1
GGACT
ADGRB2
LINGO2
TMEM109
RMDN3
RNF165
CHST2
LCNL1
NCALD
SH2D1B
TTC38
TRGV9
NEIL1
LGR6
AKR1C3
MVD
S1PR5
KLRB1
KLRD1
COPZ2
ADGRG1
BEND3
CTU2
ENSG00000261888
YPEL1
PTGDR
FCRL6
HOPX
PRSS23
KANK3
TPST2
LAMA2
PDGFD
PLEKHF1
ABCA2
ACSL5
ITPRIPL1
PAFAH2
PTGDS
PCDH1
MATK
TBX21
FGFBP2
ELOVL6
IL18RAP
BNC2
GPR68
AXIN1
DHX35
GNPTAB
CCDC28B
NCR3
B3GAT1
LLGL2
GPR153
ST8SIA6
PDGFRB
NCR1
GNRHR2
PRR5
TMEM241
ATXN7L3B
AGAP1
NKG7
DLG5
AGK
HENMT1
C1orf21
PLEKHG3
FAM131B
GZMM
TOX
GLOD4
SCD5
PPP2R2B
ENPP5
YES1
KLRK1
HIRA
GZMA
B4GALT6
SYT11
PPM1L
LINC02084
RUNX3
DTHD1
SGSM1
C19orf12
ENSG00000245025
ENPP4
VANGL1
NELFCD
PCCB
RFTN1
ANXA6
CST7
PRSS30P
CACNA2D2
GZMB
ARNTL
GNGT2
LGALS9C
PTPRA
APMAP
RGS3
SLAMF7
MIR9-3HG
MYBL1
CFH
RAP1GAP2
ABCB1
NCR3LG1
FAT4
KLRC4-KLRK1
TGFBR3
RHEBL1
CMKLR1
SLC9A3R1
ADHFE1
FUT11
PTPN4
COLQ
ENSG00000261754
FASLG
ALDH1B1
TBC1D31
TNFRSF14-AS1
TRAPPC10
SLC24A1
MBP
TSEN54
BMERB1
LINC01801
ENSG00000286010
ATP2B4
XCL2
HMOX2
JAK1
GSE1
STOM
CTSW
PRR5L
TRPV3
HLCS
APOBEC3G
IRAK2
ADGRG5
FKBP11
AKAP5
RTKN
KLRC4
JAKMIP2
DRAXIN
CRY1
GK5
TTLL1
TOGARAM2
RASSF1
LPCAT1
NOP14-AS1
GALNT11
SAMD3
IL18R1
PTPN22
FYN
PIP4K2A
PYHIN1
PHLDB2
PIP4K2C
RAB11FIP4
RAB37
TFDP2
SLC4A4
POLR3A
CD160
PPP2R5C
MAPRE2
SH2D2A
COLGALT2
CEP78
PTCH1
SPRING1
ARIH2
PIK3R3
SSH1
TSHZ1
DCLRE1A
YARS1
CHST12
STK39
TIGIT
EOMES
RHOBTB3
CAPN12
ADRB2
ENSG00000167912
NPC1
CYTH3
DERL1
LINC01145
ABHD17C
ZBTB42
PLAAT4
B3GLCT
GOLGA8N
IL12RB1
F2R
ASCL2
FYCO1
RAG1
USP28
PAXX
SYNE1
HIPK2
RNF115
ENSG00000258875
SYTL3
IL12RB2
TLR3
MTMR2
CCDC15
EPB41L4A
SYTL2
CMC1
CENPO
CCNJL
ERMP1
PCNX4
ABHD15
ZNF600
ZC3H12A-DT
LRRC1
PRIMPOL"



Cyt_genes <- strsplit(gene_list, "\\s+")[[1]]
Cyt_genes


# Filter rows where gene_name is in selected_genes
sig_genes <- df[df$gene_name %in% Cyt_genes, ]
sig_genes

sig_genes<- sig_genes[,-1]

sig_genes


#### flip this to merge it with the meta
# Step 1: Set gene names as rownames
rownames(sig_genes) <- sig_genes$gene_name

# Step 2: Drop the 'gene_name' column to keep only numeric values
sig_genes <- sig_genes[, -which(names(sig_genes) == "gene_name")]

# Step 3: Transpose the data
sig_genes_t <- as.data.frame(t(sig_genes))

sig_genes_t
# Step 4: Optional - make gene names proper column names
#colnames(df_t) <- rownames(df)

# Step 5: Add sample names as a column if you want
sig_genes_t$anon_sample_ID <- rownames(sig_genes_t)



cyto_module_genes_meta<- merge(sig_genes_t, inflm_meta, by= "anon_sample_ID")


cyto_module_genes_meta

cyto_module_genes_meta

colnames(cyto_module_genes_meta)

##### want to find the average for each sample 
df_summary <- cyto_module_genes_meta %>%
  rowwise() %>%
  dplyr::mutate(
    Cytotoxic_mean = mean(c_across(CFH:LINC01801), na.rm = TRUE),
    Cytotoxic_median = median(c_across(CFH:LINC01801), na.rm = TRUE)
  ) %>%
  ungroup()


df_summary

colnames(df_summary)

df_summary$Time_inf<- paste(df_summary$Timepoint, df_summary$inflam.gs, sep = "_")

df_summary$Time_inf<- factor(df_summary$Time_inf, levels = c( "T0_low", "T1_low", "T0_high", "T1_high")) 

df_summary <- df_summary %>%
  filter(!is.na(inflam.gs))


##### want to find the average for each sample 

my_comparisons_2 <- list(
  #c("T0_low", "T1_low"),
  c("T0_low", "T1_low"),
  c("T0_high", "T1_high"),
  c("T0_low", "T0_high"),
  c("T1_low", "T1_high"))


table(df_summary$Cohort_time)

df_summary_PCV<- df_summary[df_summary$Treatment == "PCV13", ]

colnames(df_summary_PCV)

colnames(coldata)[which(names(coldata) == "Cohort_time")] <- "Timepoint by Treatment"


flow_plot<-ggplot(df_summary_PCV, aes(x= Time_inf, y= `Cytotoxic_median`)) +
  geom_boxplot(aes(fill = Time_inf), width = 0.35,
               outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill = Time_inf), shape = 21,
              width = 0.25, alpha = 1, size = 3) +
  scale_fill_manual(
    name = "Inflammation tier",
    values = c(
      "T0_low"  = "#FFEB3B",
      "T1_low"  = "#FFEB3B",
      "T0_high" = "red",
      "T1_high" = "red"
    ),
    labels = c(
      "T0_low"  = "T0 (Low)",
      "T1_low"  = "T1 (Low)",
      "T0_high" = "T0 (High)",
      "T1_high" = "T1 (High)"
    )
  ) +
  scale_x_discrete(labels = c(
    "T0_low"  = "T0 (Low)",
    "T1_low"  = "T1 (Low)",
    "T0_high" = "T0 (High)",
    "T1_high" = "T1 (High)"
  )) +
  theme_classic() +
  ylab("Module Score")+
  xlab("")+
  ggtitle("Cytoxic gene module")+
  theme(plot.title = element_text(hjust=0.5))+
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black")  +
  #geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comparisons_2, method= "t.test",  label = "p.signif", hide.ns = FALSE )+
  theme(
    #legend.position = "none", 
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  )
#stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comparisons_2, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE )# + #label.y = c(95, 90, 85, 80, 75, 40) 


flow_plot

ggsave(plot = flow_plot, "Figures/Modules/Paper/Cytotoxic_gene_score.png", dpi = 1200,  height = 9,  width = 9)


####### plasmablast module #######
gene_list_2<-"TXNDC5 
TNFRSF17
LOC652694
LOC647450
LOC651751
TXNDC5
LOC652493
LOC642113
GLDC
LOC652102
LOC649210
LOC649923
CD38
JCHAIN
LOC652113
LOC647460
LOC652775
LOC401845
CLICKIII"



plas_genes <- strsplit(gene_list_2, "\\s+")[[1]]
plas_genes

df
# Filter rows where gene_name is in selected_genes
sig_genes_2 <- df[df$gene_name %in% plas_genes, ]
sig_genes_2

sig_genes_2<- sig_genes_2[,-1]

sig_genes_2


#### flip this to merge it with the meta
# Step 1: Set gene names as rownames
rownames(sig_genes_2) <- sig_genes_2$gene_name

# Step 2: Drop the 'gene_name' column to keep only numeric values
sig_genes_2 <- sig_genes_2[, -which(names(sig_genes_2) == "gene_name")]

# Step 3: Transpose the data
sig_genes_2_t <- as.data.frame(t(sig_genes_2))

sig_genes_2_t
# Step 4: Optional - make gene names proper column names
#colnames(df_t) <- rownames(df)

# Step 5: Add sample names as a column if you want
sig_genes_2_t$anon_sample_ID <- rownames(sig_genes_2_t)

sig_genes_2_t

plasma_module_genes_meta<- merge(sig_genes_2_t, inflm_meta, by= "anon_sample_ID")


plasma_module_genes_meta


colnames(plasma_module_genes_meta)

##### want to find the average for each sample 
df_summary_2 <- plasma_module_genes_meta %>%
  rowwise() %>%
  dplyr::mutate(
    Plasmablast_mean = mean(c_across(CD38:TXNDC5), na.rm = TRUE),
    Plasmablast_median = median(c_across(CD38:TXNDC5), na.rm = TRUE)
  ) %>%
  ungroup()


df_summary_2

colnames(df_summary_2)

df_summary_2$Time_inf<- paste(df_summary_2$Timepoint, df_summary_2$inflam.gs, sep = "_")

df_summary_2

df_summary_2<- df_summary_2[df_summary_2$Treatment == "PCV13",]

df_summary_2
colnames(df_summary_2)

df_summary_2_clean <- df_summary_2 %>% filter(!is.na(inflam.gs))
df_summary_2_clean$Time_inf<- factor(df_summary_2_clean$Time_inf, levels = c( "T0_low", "T1_low", "T0_high", "T1_high")) 

my_comparisons_2 <- list(
  #c("T0_low", "T1_low"),
  c("T0_low", "T1_low"),
  c("T0_high", "T1_high"),
  c("T0_low", "T0_high"),
  c("T1_low", "T1_high"))



flow_plot<-ggplot(df_summary_2_clean, aes(x= Time_inf, y= `Plasmablast_median`)) +
  geom_boxplot(aes(fill = Time_inf), width = 0.35,
               outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(fill = Time_inf), shape = 21,
              width = 0.25, alpha = 1, size = 3) +
  scale_fill_manual(
    name = "Inflammation tier",
    values = c(
      "T0_low"  = "#FFEB3B",
      "T1_low"  = "#FFEB3B",
      "T0_high" = "red",
      "T1_high" = "red"
    ),
    labels = c(
      "T0_low"  = "T0 (Low)",
      "T1_low"  = "T1 (Low)",
      "T0_high" = "T0 (High)",
      "T1_high" = "T1 (High)"
    )
  ) +
  scale_x_discrete(labels = c(
    "T0_low"  = "T0 (Low)",
    "T1_low"  = "T1 (Low)",
    "T0_high" = "T0 (High)",
    "T1_high" = "T1 (High)"
  )) +
  theme_classic() +
  ylab("Module Score")+
  xlab("")+
  ggtitle("Plasmablast gene module (M4.11)")+
  theme(plot.title = element_text(hjust=0.5))+
  geom_vline(xintercept = 2.5, linetype= "dotted", color = "black")  +
  #geom_vline(xintercept = 4.5, linetype= "dotted", color = "black")  +
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comparisons_2, method= "t.test",  label = "p.signif", hide.ns = FALSE )+
  theme(
    #legend.position = "none", 
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1.2, "cm")
  )
#stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comparisons_2, method= "wilcox.test",  label = "p.signif", hide.ns = FALSE )# + #label.y = c(95, 90, 85, 80, 75, 40) 


flow_plot

ggsave(plot = flow_plot, "Figures/Modules/Paper/Plasmablast_gene_score.png", dpi = 1200,  height = 9,  width = 9)



######       Volcano Plot       ########
#####
####
###
##
#

#load df
df <- read.csv("Transcript/new_analysis/results_T0_PCVInflamhivinflamlow.csv")
df <- read.csv("Transcript/new_analysis/results_T0_PCVvPLacebo_withInflamtiers.csv")

## add a column of NAs
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$log2FoldChange > 1 & df$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -1 & df$padj < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$gene_name[df$diffexpressed != "NO"]


mycolors <- c("dodgerblue", "firebrick1", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")


####### filter the df to only include genes from the cytotoxic module 


colnames(df)
df
# Filter rows where gene_name is in selected_genes
cyto_genes <- df[df$gene_name %in% Cyt_genes, ]



VP_df <- ggplot(data=cyto_genes, aes(x=(log2FoldChange), y=-log10(padj), label=delabel))+ 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 5)+
  theme_classic()+
  geom_vline(xintercept=c(-1, 1), linetype="dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", col="black")+
  scale_fill_manual(values=mycolors)+
  geom_label_repel(size = 10, force = 2, max.overlaps = 10)+
  ggtitle("PCV13 T0 v Placebo T0")+
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
  theme(plot.title = element_text(hjust=0.5, face = "bold"))+
  ylim(0,15)+
  xlim(-3, 3)

VP_df
ggsave(plot = VP_df, "Figures/Modules/paper/VP_inflamhivinflamlo_PCVT0vPlacebo_cytotoxicgenes.png", dpi = 1200, height = 9, width = 9, units = "in")

###########################

inflm_meta<- read.csv("Gitools/df_annotations_2tiers.csv")
inflm_meta



#### bring in serology results 

WHO_T0_T2<- read.csv("paper/final anonymised dataframes/Serology/serology_anon.csv")
WHO_T0_T2

#### change above back to the WHO df
colnames(WHO_T0_T2)
colnames(inflm_meta)

colnames(inflm_meta)[which(names(inflm_meta) == "X")] <- "anon_sample_ID"


inflm_meta_sero<- merge( WHO_T0_T2 ,inflm_meta, by="anon_sample_ID", all.x=TRUE)

inflm_meta_sero


inflm_meta_sero$inflam.gs

inflm_meta_sero_baseline <- inflm_meta_sero %>% filter(!is.na(inflam.gs))

inflm_meta_sero_baseline
colnames(inflm_meta_sero_baseline)
inflm_meta_sero_baseline$Clusters<- as.factor(inflm_meta_sero_baseline$Clusters)

my_comps_3<- list( c("tier1", "tier2"), c("tier1", "tier3"), c("tier2", "tier3"))
my_comps_3<- list( c("1", "2"))

dev
################### boxplot Serotpye with cohort_time ######

colnames(inflm_meta_sero_baseline)

inflm_meta_sero$timepoint

serotype<-ggplot(inflm_meta_sero_baseline, aes(x= inflam.gs, y= IgG_type.23F)) +
  geom_boxplot(aes(fill= inflam.gs), width=0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= inflam.gs), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("low" = "yellow",
                              #"tier2" ="orange",
                              "high" = "red"))+
  #scale_fill_manual(values= c("1" = "blue",
  #                            "2" ="green"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("IgG ug/ml")+
  xlab("")+
  ggtitle("23F")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  #geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = list( c("PCV13 T0", "PCV13 T2")), method= "annova",  label = "p.signif", hide.ns = TRUE)
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_3, label = "p.signif", hide.ns = FALSE, method = "wilcox.test")
#geom_signif(comparisons = my_comps_3, map_signif_level = TRUE, test = "t.test")

serotype
ggsave(plot = serotype, "Figures/Gitools/Baseline_23F.png", dpi = 300,  height = 6,  width = 6)

#####


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
  "Figures/Gitools/inflam_legend.png",
  legend_only,
  width = 10,
  height = 2,
  dpi = 1200
)

#c("Placebo T0", "Placebo T2"),
#ggsave(plot = serotype, "Figures/Serology/Serotype23F_box.png", dpi = 300,  height = 6,  width = 6)
########

inflm_meta_sero



####################

inflm_meta_sero


inflm_meta_sero_PCV<- inflm_meta_sero[inflm_meta_sero$group == "PCV13", ]
inflm_meta_sero_PCV
colnames(inflm_meta_sero_PCV)

# Step 1: Extract inflam.gs values from T0
inflam_T0 <- inflm_meta_sero_PCV %>%
  filter(timepoint == "T0") %>%
  select(patient_num, inflam.gs_T0 = inflam.gs, Cluster_T0 = Clusters) 

# Step 2: Left join onto the full dataframe by patient_num
inflm_meta_sero_PCV_full <- inflm_meta_sero_PCV %>%
  left_join(inflam_T0, by = "patient_num") %>%
  mutate(inflam.gs = ifelse(is.na(inflam.gs), inflam.gs_T0, inflam.gs)) %>%
  mutate(Clusters = ifelse(is.na(Clusters), Cluster_T0, Clusters)) %>%
  select(-inflam.gs_T0, -Cluster_T0)  

inflm_meta_sero_PCV_full

inflm_meta_sero_PCV_full$Time_inflam<- paste(inflm_meta_sero_PCV_full$timepoint, inflm_meta_sero_PCV_full$inflam.gs, sep = "_")
inflm_meta_sero_PCV_full$Time_cluster<- paste(inflm_meta_sero_PCV_full$timepoint, inflm_meta_sero_PCV_full$Clusters, sep = "_")

inflm_meta_sero_PCV_full

table(inflm_meta_sero_PCV_full$Time_inflam)

inflm_meta_sero_PCV_full <- inflm_meta_sero_PCV_full %>% filter(!is.na(inflam.gs))
inflm_meta_sero_PCV_full


#inflm_meta_sero_PCV_full$Time_inflam<- factor(inflm_meta_sero_PCV_full$Time_inflam, levels = c("T0_tier1", "T2_tier1", "T0_tier2","T2_tier2","T0_tier3" ,"T2_tier3")) 

inflm_meta_sero_PCV_full$Time_inflam<- factor(inflm_meta_sero_PCV_full$Time_inflam, levels = c("T0_low", "T2_low", "T0_high","T2_high")) 


my_comps_3<- list( c("T0_tier1", "T2_tier1"),c("T0_tier2", "T2_tier2"),c("T0_tier3", "T2_tier3"))


colnames(inflm_meta_sero_PCV_full)


# Split data by species
iris_split <- split(inflm_meta_sero_PCV_full$IgG_type.23F, inflm_meta_sero_PCV_full$Time_inflam)
iris_split
# Apply Shapiro-Wilk test to each group
lapply(iris_split, shapiro.test) ### if greater than 0.05 then it is not normal 



############## to do them individually 
t_test_results <- compare_means(IgG_type.6B ~ Time_inflam, data = inflm_meta_sero_PCV_full, method = "wilcox.test")
t_test_results


significant_results <- t_test_results %>% 
  filter(p < 0.05)  # Adjust the threshold if needed

significant_results


serotype<-ggplot(inflm_meta_sero_PCV_full, aes(x= Time_inflam, y= IgG_type.23F)) +
  geom_boxplot(aes(fill= Time_inflam), width=0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Time_inflam), width=0.25, shape = 21,alpha=1, size = 4) +
  scale_fill_manual(values= c("T0_low" = "#fbd506",
                              "T2_low" = "#fbd506",
                              #"T0_tier2" ="#FA9600",
                              #"T2_tier2" ="#FA9600",
                              "T0_high" = "red", 
                              "T2_high" = "red"))+
  #scale_fill_manual(values= c("1" = "blue",
  #                            "2" ="green"))+
  theme_classic()+
  scale_y_log10()+ # limits= c(0.1, 10)
  annotation_logticks(sides = "l")+
  ylab("IgG ug/ml")+
  xlab("")+
  ggtitle("23F")+
  theme(plot.title = element_text(hjust=0.5))+
  #geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
  #geom_vline(xintercept = 2.5, linetype= "dotted", color = "black") #+
  #stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = list( c("PCV13 T0", "PCV13 T2")), method= "annova",  label = "p.signif", hide.ns = TRUE)
  stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_3, label = "p.signif", hide.ns = FALSE, method = "wilcox.test") +
  scale_x_discrete(labels=c('T0 low', 'T2 low', 'T0 high', 'T2 high'))
#geom_signif(comparisons = my_comps_3, map_signif_level = TRUE, test = "t.test")

serotype

ggsave(plot = serotype, "Figures/Gitools/inflm_T0T2_2tiers_23F.png", dpi = 300,  height = 6,  width = 6)





