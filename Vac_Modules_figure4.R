######### Figure 3 Modules #######

library(tidyverse)
library(tximport)
library(readr)
library(GenomicFeatures)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(VennDiagram)
#library(ggvenn)
#library(eulerr)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(pheatmap)
library(ComplexHeatmap)
library(viridis)
#BiocManager::install("qusage", force = TRUE)


#########


##### bring in the normalised output of deseq 

counts<- read.csv("paper/final anonymised dataframes/Transcript/Normalisedcounts_ddsDeSEQ2_Cohort_time_anon.csv") 
counts

counts<- counts %>% column_to_rownames("gene_name")
rownames(counts)<- make.names(counts$gene_name, unique = TRUE)

counts<- counts[,-c(1, 2)]
### might need to add 1 before logging to get rid of infinites 

counts<- counts +1

counts_log<- log2(counts)
counts
counts_log

#load metadata
metadata <- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata.csv") #400
metadata

colnames(metadata)
colnames(counts_log)
sample_names <- colnames(counts_log)
sample_names
duplicated(sample_names)

# Filter metadata to only keep rows where sample_ID is in the expression data
meta_filtered <- metadata %>%
  dplyr::filter(anon_sample_ID %in% sample_names)

metadata<- meta_filtered %>% dplyr::select(patient_num, anon_sample_ID, Timepoint, Treatment, Cohort_time)

metadata  


metadata<- metadata %>% column_to_rownames("anon_sample_ID")


all(rownames(metadata) %in% colnames(counts_log))
setdiff(rownames(metadata), colnames(counts_log))

metadata$anon_sample_ID<- rownames(metadata)


# Reorder columns to match a specific ID column
counts_log <- counts_log[, match(metadata$anon_sample_ID, colnames(counts_log))]


metadata

######### labels


labels<- metadata$Cohort_time
labels

#### pairs


pairs<- metadata$patient_num
pairs

##### gene sets 
MSIG.geneSets = read.gmt("Transcript/h.all.v2024.1.Hs.symbols.gmt")
MSIG.geneSets$HALLMARK_INFLAMMATORY_RESPONSE

##### gene sets 
BTM.geneSets = read.gmt("Transcript/BTM_for_GSEA_20131008.gmt")
BTM.geneSets

Reactome.geneSets = read.gmt("Transcript/c2.cp.reactome.v2024.1.Hs.symbols.gmt")

Reactome.geneSets$REACTOME_INFLAMMASOMES


Epithelial.geneSets = read.gmt("Transcript/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt")
Epithelial.geneSets

summary(MSIG.geneSets[1:20])

summary(BTM.geneSets[1:20])

####### Qusage ##############

############ contrast PCV timepoints


contrast = "PCV13_T1-PCV13_T0"


contrast
qs.results_PCV = qusage(counts_log, labels, contrast, BTM.geneSets, pairVector = pairs)
###### warning message GeneSet 'interferon alpha response (II) (M158.1)' contains one or zero overlapping genes. NAs produced. ##### probably best not to use it? 

### 
p.vals = pdf.pVal(qs.results_PCV)


head(p.vals)
is.vector(p.vals)

q.vals = p.adjust(p.vals, method="fdr")

head(q.vals)
class(q.vals)

table_PCV_T1<-qsTable(qs.results_PCV, number=1000, sort.by = c("fdr"))
table_PCV_T1
summary(table_PCV_T1$log.fold.change)


class(table_PCV_T1)
class(table_PCV_T1_hall)



write.csv(table_PCV_T1, "Transcript/Qusage/T1PCV_resultstable.csv")



plot(qs.results_PCV, xlab="Gene Set Activation")

pdf("Figures/Qusage/gene_set_activation_plot_T1PCV.pdf", width = 20, height = 20)  # Width and height in inches
plot(qs.results_PCV, xlab = "PCV13 T1-T0 Gene Set Activation")
dev.off()

numPathways(qs.results_PCV)

qs.results_PCV
#### mean activity scores 
qs.results_PCV$
  mean_scores_PCVT1_base<-qs.results_PCV$path.mean
mean_scores_PCVT1_base
class(mean_scores_PCVT1_base)
mean_scores_PCVT1_base

mean_scores_PCVT1_base_df<- as.data.frame(mean_scores_PCVT1_base)
mean_scores_PCVT1_base_df

write.csv(mean_scores_PCVT1_base_df, "Transcript/Qusage/PCV13T1-T0_meanactivityscores.csv")


###### PCV13 T2- PCV13 T0 #################

contrast = "PCV13_T2-PCV13_T0"

qs.results_PCVT2 = qusage(counts_log, labels, contrast, BTM.geneSets, pairVector = pairs)


p.vals_T2 = pdf.pVal(qs.results_PCVT2)
head(p.vals_T2)


q.vals_T2 = p.adjust(p.vals_T2, method="fdr")
head(q.vals_T2)
q.vals_T2

table_PCV_T2<-qsTable(qs.results_PCVT2, number=1000, sort.by = c("fdr"))
table_PCV_T2


write.csv(table_PCV_T2, "Transcript/Qusage/T2PCV_resultstable.csv")




pdf("Figures/Qusage/gene_set_activation_plot_T2PCV.pdf", width = 20, height = 20)  # Width and height in inches
plot(qs.results_PCVT2, xlab = "PCV13 T2 - PCV13 T0 Gene Set Activation")
dev.off()

#### mean activity scores 

mean_scores_T2PCV<-qs.results_PCVT2$path.mean
class(mean_scores_T2PCV)


mean_scores_T2PCV_df<- as.data.frame(mean_scores_T2PCV)
mean_scores_T2PCV_df

write.csv(mean_scores_T2PCV_df, "Transcript/Qusage/PCV13T2-T0_meanactivityscores.csv")




#####################
###### Placebo ###### 
######         ###### 


contrast = "Placebo_T1-Placebo_T0"

qs.results_placebo = qusage(counts_log, labels, contrast, BTM.geneSets, pairVector = pairs)

p.vals_placebo = pdf.pVal(qs.results_placebo)
head(p.vals_placebo)


q.vals_placebo = p.adjust(p.vals_placebo, method="fdr")
head(q.vals_placebo)
q.vals_placebo

table_placebo<-qsTable(qs.results_placebo, number=1000, sort.by = c("fdr"))
table_placebo

write.csv(table_placebo, "Transcript/Qusage/T1-Placebo_resultstable.csv")



pdf("Figures/Qusage/gene_set_activation_plot_T1Placebo.pdf", width = 20, height = 20)  # Width and height in inches
plot(qs.results_placebo, xlab = "Placebo T1 - Placebo T0 Gene Set Activation")
dev.off()


class(q.vals_placebo)

class(table)
table_placebo
write.csv(table_placebo, "Transcript/Qusage/T1-Placebo_resultstable.csv")

#### mean activity scores 

mean_scores_PlaceboT1_base<-qs.results_placebo$path.mean
class(mean_scores_PlaceboT1_base)
mean_scores_PlaceboT1_base

mean_scores_PlaceboT1_base_df<- as.data.frame(mean_scores_PlaceboT1_base)
mean_scores_PlaceboT1_base_df

write.csv(mean_scores_PlaceboT1_base_df, "Transcript/Qusage/PlaceboT1-T0_meanactivityscores.csv")


################# ################# #################
###### Placebo T2- Placebo T0 #################
################# ################# #################
contrast = "Placebo_T2-Placebo_T0"

qs.results_PlaceboT2 = qusage(counts_log, labels, contrast, BTM.geneSets, pairVector = pairs)

p.vals_T2 = pdf.pVal(qs.results_PlaceboT2)
head(p.vals_T2)


q.vals_T2 = p.adjust(p.vals_T2, method="fdr")
head(q.vals_T2)
q.vals_T2

table_T2<-qsTable(qs.results_PlaceboT2, number=1000, sort.by = c("fdr"))
table_T2


write.csv(table_T2, "Transcript/Qusage/T2Placebo_resultstable.csv")



pdf("Figures/Qusage/gene_set_activation_plot_T2Placebo.pdf", width = 20, height = 20)  # Width and height in inches
plot(qs.results_PlaceboT2, xlab = "Placebo T2 - Placebo T0 Gene Set Activation")
dev.off()

#### mean activity scores 

mean_scores_T2Placebo<-qs.results_PlaceboT2$path.mean
class(mean_scores_T2Placebo)


mean_scores_T2Placebo_df<- as.data.frame(mean_scores_T2Placebo)
mean_scores_T2Placebo_df

write.csv(mean_scores_T2Placebo_df, "Transcript/Qusage/PClaceboT2-T0_meanactivityscores.csv")

##################################################################
########### heatmap prep  ###########

### bring in the results tables ####

Placebo_1<- read.csv("Transcript/Qusage/T1-Placebo_resultstable.csv")
Placebo_2<- read.csv("Transcript/Qusage/T2Placebo_resultstable.csv")


PCV_1<- read.csv("Transcript/Qusage/T1PCV_resultstable.csv")
PCV_2<- read.csv("Transcript/Qusage/T2PCV_resultstable.csv")

#### now merge all the original unfiltered tables together 

colnames(Placebo_1)[which(names(Placebo_1) == "log.fold.change")] <- "Placebo_T1"
colnames(Placebo_1)[which(names(Placebo_1) == "FDR")] <- "FDR_Placebo_T1"
Placebo_1<- Placebo_1[,-c(1,4)]
Placebo_1

colnames(Placebo_2)[which(names(Placebo_2) == "log.fold.change")] <- "Placebo_T2"
colnames(Placebo_2)[which(names(Placebo_2) == "FDR")] <- "FDR_Placebo_T2"
Placebo_2<- Placebo_2[,-c(1,4)]
Placebo_2



colnames(PCV_1)[which(names(PCV_1) == "log.fold.change")] <- "PCV13_T1"
colnames(PCV_1)[which(names(PCV_1) == "FDR")] <- "FDR_PCV13_T1"
PCV_1<- PCV_1[,-c(1,4)]
PCV_1


colnames(PCV_2)[which(names(PCV_2) == "log.fold.change")] <- "PCV13_T2"
colnames(PCV_2)[which(names(PCV_2) == "FDR")] <- "FDR_PCV13_T2"
PCV_2<- PCV_2[,-c(1,4)]
PCV_2


###bind them together 

# List of original results tables dataframes to merge
dfs <- list(PCV_1, PCV_2, Placebo_1, Placebo_2 )

PCV_1

merged_df <- Reduce(function(x, y) merge(x, y, by = "pathway.name", all = TRUE), dfs)

merged_df
write.csv(merged_df, "Transcript/Qusage/Merged_resultstables_Tx-T0_anon.csv")
##########
merged_df<- read.csv("Transcript/Qusage/Merged_resultstables_Tx-T0_anon.csv")

merged_df <- merged_df[,-1]

merged_df


############ make a df of just B cell associated pathways 
### use the big df not the signifpathways one


merged_df_scores<- merged_df  %>% dplyr::select(-starts_with("FDR_"))
merged_df_scores
merged_df_scores_m <- merged_df_scores %>% column_to_rownames("pathway.name")
merged_df_scores_m<- as.matrix(merged_df_scores_m)
colnames(merged_df_scores_m)

merged_df_scores_m <- merged_df_scores_m[, c("PCV13_T1", "Placebo_T1", "PCV13_T2", "Placebo_T2" )]


matches <- BTM.geneSets[grep("B cell.*|Plasma.*|Plasmablast.*|Antibody.*|Immunoglobulins.*|BCR.*", names(BTM.geneSets), ignore.case = TRUE)]

matches

b_vector<- names(matches)
b_vector

Bcell_mat <- merged_df_scores_m[rownames(merged_df_scores_m) %in% b_vector, ]

Bcell_mat


# Rowname to remove
remove_row <- c("plasma membrane, cell junction (M162.0)","signal transduction, plasma membrane (M82)", "enriched in plasma membrane proteins (II) (M135.1)", "enriched in plasma membrane proteins (I) (M135.0)", "immuregulation - monocytes, T and B cells (M57)")


# Subset the matrix, excluding rows with rownames in `remove_row`
Bcell_mat_trim <- Bcell_mat[!rownames(Bcell_mat) %in% remove_row, ]



pheatmap(Bcell_mat_trim, scale = "none", cluster_cols =F , cluster_rows = T)


###### try and get a significance column 

colnames(merged_df)
merged_df
merged_df_T0_T2<- merged_df

merged_df_sig<- merged_df_T0_T2

merged_df_sig
colnames(merged_df_sig)

# Define the significance threshold
threshold <- 0.05

# Create a matrix of stars based on the FDR values
star_matrix <- matrix("", nrow = nrow(merged_df_sig), ncol = 4)  # Adjust the number of columns as per heatmap
rownames(star_matrix) <- merged_df_sig$pathway.name
colnames(star_matrix) <- c("PCV13_T1", "PCV13_T2",  "Placebo_T1", "Placebo_T2")


# Fill in stars for significant FDR values
star_matrix[, "PCV13_T1"] <- ifelse(merged_df_sig$FDR_PCV13_T1 < threshold, "*", "")
star_matrix[, "PCV13_T2"] <- ifelse(merged_df_sig$FDR_PCV13_T2 < threshold, "*", "")
star_matrix[, "Placebo_T1"] <- ifelse(merged_df_sig$FDR_Placebo_T1 < threshold, "*", "")
star_matrix[, "Placebo_T2"] <- ifelse(merged_df_sig$FDR_Placebo_T2 < threshold, "*", "")


star_matrix

Bcell_mat_trim

Bcell_mat_trim_T0_T2<- Bcell_mat_trim

Bcell_mat_trim_T0_T2

dim(Bcell_mat_trim_T0_T2)
dim(star_matrix)

b_vector

Bcell_mat_trim_T0_T2


star_matrix_B <- star_matrix[rownames(star_matrix) %in% b_vector, ]

star_matrix_B


star_matrix_B <- star_matrix_B[!rownames(star_matrix_B) %in% remove_row, ]

dim(Bcell_mat_trim_T0_T2)
dim(star_matrix_B)

# Plot the heatmap with pheatmap and add the stars as annotations
pheatmap(
  Bcell_mat_trim_T0_T2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = star_matrix_B,  # Overlay stars for significant FDRs
  number_color = "black",         # Color of the stars
  fontsize_number = 16,
  scale = "row",
  main = "Heatmap with Significant FDRs",
  cellwidth = 20
)
Bcell_mat_trim_T0_T2

colnames(Bcell_mat_trim_T0_T2) <- c("PCV13 T1 v T0", "Placebo T1 v T0", "PCV13 T2 v T0", "Placebo T2 v T0")
library(grid)


pdf("Figures/Qusage/Heatmap_B_modules_significance_0.05_anon.pdf", height = 20 , width = 20)
png("Figures/Qusage/Heatmap_B_modules_significance_0.05_anon_legend.png", height = 15 , width = 15,  units = "in", res = 700)

pheatmap::pheatmap(
  Bcell_mat_trim_T0_T2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = star_matrix_B,  # Overlay stars for significant FDRs
  number_color = "black",         # Color of the stars
  fontsize = 25,
  fontsize_number = 20,
  fontsize_col = 26,
  fontsize_row = 28,
  scale = "none",
  main = "", #B cell modules
  cellwidth = 25,
  cellheight = 35,
  legend = TRUE,
  height = 20,
  width = 20
)
dev.off()

library(ComplexHeatmap)
library(grid) # for gpar
library(circlize)

# Your custom palette
min_val <- min(Bcell_mat_trim_T0_T2, na.rm = TRUE)
min_val
max_val <- max(Bcell_mat_trim_T0_T2, na.rm = TRUE)
max_val
my_col_fun <- colorRamp2(
  seq(min_val, max_val, length.out = 5),
  c("#163f76", "#4575B4", "#FFFFF2", "#ec763f", "#bf2a1b")
)


Heatmap(
  Bcell_mat_trim_T0_T2,
  name = "B cell modules",
  col = my_col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_matrix_B[i, j], x, y,
              gp = gpar(fontsize = 16, col = "black"))
  },
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_title = "B cell modules"
)



######################################
### T cells ##### 
#######################################

matches_T <- BTM.geneSets[grep("T cell.*|CD4.*|CD8.*|TCR.*", names(BTM.geneSets), ignore.case = TRUE)]

matches_T

T_vector<- names(matches_T)
T_vector

Tcell_mat <- merged_df_scores_m[rownames(merged_df_scores_m) %in% T_vector, ]

Tcell_mat


pheatmap(Tcell_mat, scale = "row", cluster_cols =T , cluster_rows = T)


#pdf("Figures/Qusage/Heatmap_timepoints_Tcells.pdf", height = 20 , width = 20)
pheatmap(Tcell_mat, scale = "none", cluster_cols =F , cluster_rows = T)
#dev.off()



######## with significance ### 

##

merged_df_sig ### this is everypathway with FDR 
star_matrix ### every pathway with * when 


### take the T cell matrix and get rid of T3
colnames(Tcell_mat)
Tcell_mat
Tcell_mat_T0_T2 <- Tcell_mat

### then need to suset the star matrix to only match the T cell one, use dim() to check 


T_vector

star_matrix_T <- star_matrix[rownames(star_matrix) %in% T_vector, ]


dim(Tcell_mat_T0_T2)
dim(star_matrix_T)

# Plot the heatmap with pheatmap and add the stars as annotations
pheatmap(
  Tcell_mat_T0_T2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = star_matrix_T,  # Overlay stars for significant FDRs
  number_color = "black",         # Color of the stars
  fontsize_number = 16,
  scale = "none",
  main = "Heatmap with Significant FDRs"
)

colnames(Tcell_mat_T0_T2) <- c("PCV13 T1 v T0", "Placebo T1 v T0", "PCV13 T2 v T0", "Placebo T2 v T0")

# Your custom palette

min_val <- min(Tcell_mat_T0_T2, na.rm = TRUE)
min_val
max_val <- max(Tcell_mat_T0_T2, na.rm = TRUE)
max_val
my_col_fun <- colorRamp2(
  seq(min_val, max_val, length.out = 5),
  c("#163f76", "#4575B4", "#FFFFF2", "#ec763f", "#bf2a1b")
)


Heatmap(
  Tcell_mat_T0_T2,
  name = "T cell modules",
  col = my_col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_matrix_T[i, j], x, y,
              gp = gpar(fontsize = 16, col = "black"))
  },
  column_names_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 14, fontface = "bold"),
  column_title = "T cell modules"
)

?pheatmap
pdf("Figures/Qusage/Heatmap_T_modules_significance_0.05_anon.pdf", height = 20 , width = 20)
png("Figures/Qusage/Heatmap_T_modules_significance_0.05_anon_legend.png", height = 15 , width = 15,  units = "in", res = 700)

pheatmap::pheatmap(
  Tcell_mat_T0_T2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = star_matrix_T,  # Overlay stars for significant FDRs
  number_color = "black",         # Color of the stars
  fontsize = 25,
  fontsize_number = 20,
  fontsize_col = 26,
  fontsize_row = 28,
  scale = "none",
  main = "", #T cell modules
  cellwidth = 25,
  cellheight = 35,
  legend = TRUE,
  height = 20,
  width = 20
)
dev.off()




###### myeloid cells 



matches_myl <- BTM.geneSets[grep(".*dendritic cell.*|.*Dendritic.*|.*myeloid.*|.*macrophage.*|.*NK cell.*|.*monocyte.*|.*neutrophil.*|.*inflammation.*|.*inflammatory.*|.*complement.*|.*DC.*|Recruitment.*|Antigen.*", names(BTM.geneSets), ignore.case = TRUE)]

matches_myl


myeloid_vector<- names(matches_myl)
myeloid_vector

myeloid_cell_mat <- merged_df_scores_m[rownames(merged_df_scores_m) %in% myeloid_vector, ]

myeloid_cell_mat



######### myeloid with significance 

colnames(myeloid_cell_mat)
myeloid_cell_mat_T0_T2 <- myeloid_cell_mat

### then need to suset the star matrix to only match the T cell one, use dim() to check 


myeloid_vector

star_matrix_myeloid <- star_matrix[rownames(star_matrix) %in% myeloid_vector, ]




dim(myeloid_cell_mat_T0_T2)
dim(star_matrix_myeloid)

# Plot the heatmap with pheatmap and add the stars as annotations
pheatmap(
  myeloid_cell_mat_T0_T2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = star_matrix_myeloid,  # Overlay stars for significant FDRs
  number_color = "black",         # Color of the stars
  fontsize_number = 16,
  scale = "none",
  main = "Heatmap with Significant FDRs"
)

colnames(myeloid_cell_mat_T0_T2) <- c("PCV13 T1 v T0", "Placebo T1 v T0", "PCV13 T2 v T0", "Placebo T2 v T0")

rows_to_remove <- c("complement activation (I) (M112.0)", 
                    "complement activation (II) (M112.1)",
                    "complement and other receptors in DCs (M40)",
                    "extracellular matrix, complement (M140)",
                    "formyl peptide receptor mediated neutrophil response (M11.2)",
                    "inflammatory response (M33)")

# Subset the matrix to keep only rows that are NOT in the `rows_to_remove` list
myeloid_cell_mat_T0_T2 <- myeloid_cell_mat_T0_T2[!rownames(myeloid_cell_mat_T0_T2) %in% rows_to_remove, ]

star_matrix_myeloid <- star_matrix_myeloid[!rownames(star_matrix_myeloid) %in% rows_to_remove, ]


myeloid_cell_mat_T0_T2

pdf("Figures/Qusage/Heatmap_myeloid_modules_significance_0.05_legend.pdf", height = 20 , width = 20)
png("Figures/Qusage/Heatmap_myeloid_modules_significance_0.05_legend.png", height = 20 , width = 20,  units = "in", res = 700)

pheatmap::pheatmap(
  myeloid_cell_mat_T0_T2,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = star_matrix_myeloid,  # Overlay stars for significant FDRs
  number_color = "black",         # Color of the stars
  fontsize = 25,
  fontsize_number = 20,
  fontsize_col = 26,
  fontsize_row = 28,
  scale = "none",
  main = "", #Myeloid cell modules
  cellwidth = 25,
  cellheight = 28,
  legend = TRUE,
  height = 20,
  width = 20
)
dev.off()