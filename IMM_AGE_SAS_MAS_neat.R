#SAS, MAS and IMMUNE AGE volcano plots
library(limma)
library(ggplot2)
#load T0 and T1 gene expression data for PCV13 and Placebo
df_T0_T1 <- read.csv("df_T0_T1.csv", row.names = 1)
head(df_T0_T1)


################################################################################
#IMMUNE AGE

#read in IMMUNE AGE gene sets
#Genes that have a positive direction along immune age
genes_pos <- readLines("genes_pos.txt")
genes_pos
genes_neg <- readLines("genes_neg.txt")
genes_neg

#read in gene expression data for these genes
mat_genes_IMMAGE <- read.csv("new/revisions/mat_genes_IMMAGE.csv", row.names = 1)
head(mat_genes_IMMAGE)

#add in metadata
metadata <- read.csv("anon_metadata.csv")
metadata <- metadata %>% filter(Timepoint == "T0" | Timepoint == "T1")
metadata <- metadata %>% dplyr::select("anon_sample_ID", "Cohort_time")

#expression matrix: genes x samples
expr <- as.matrix(mat_genes_IMMAGE)

dim(expr)
nrow(expr)  #should be number of genes/features
ncol(expr)  #should be number of samples

#ensure sample order matches metadata
expr <- expr[, metadata$anon_sample_ID]
#samples in metadata but not in expression matrix
setdiff(metadata$anon_sample_ID, colnames(expr))
#samples in expression matrix but not in metadata
setdiff(colnames(expr), metadata$anon_sample_ID)

#Contrasts for DEGS
metadata$Cohort_time <- factor(metadata$Cohort_time)
design <- model.matrix(~ 0 + Cohort_time, data = metadata)
colnames(design) <- sub("^Cohort_time", "", colnames(design))
colnames(design)

colnames(design)
fit <- lmFit(expr, design)
contrasts <- makeContrasts(
  T0_PCV13_vs_Placebo = PCV13_T0 - Placebo_T0,
  T1_PCV13_vs_Placebo = PCV13_T1 - Placebo_T1,
  T1_Placebo_vs_T0_Placebo = Placebo_T1 - Placebo_T0,
  T1_PCV13_vs_T0_PCV13 = PCV13_T1 - PCV13_T0,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
res_T0 <- topTable(fit2, coef = "T0_PCV13_vs_Placebo", number = Inf)
res_T1 <- topTable(fit2, coef = "T1_PCV13_vs_Placebo", number = Inf)
res_Placebo <- topTable(fit2, coef = "T1_Placebo_vs_T0_Placebo", number = Inf)
res_PCV13 <- topTable(fit2, coef = "T1_PCV13_vs_T0_PCV13", number = Inf)

###########################  PCV13 T0 vs Placebo T0 ############################

res_T0 <- topTable(fit2, coef = "T0_PCV13_vs_Placebo", number = Inf) %>% rownames_to_column("Gene_name")
# add a column of NAs
res_T0$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_T0$diffexpressed[res_T0$logFC > 0.5 & res_T0$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_T0$diffexpressed[res_T0$logFC < -0.5 & res_T0$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res_T0$delabel <- NA
res_T0$delabel[res_T0$diffexpressed != "NO"] <- res_T0$Gene_name[res_T0$diffexpressed != "NO"]

res_T0$Direction_IMM_AGE <- "Other"
res_T0$Direction_IMM_AGE[res_T0$Gene_name %in% genes_pos_fixed] <- "> IMM_AGE"
res_T0$Direction_IMM_AGE[res_T0$Gene_name %in% genes_neg_fixed] <- "< IMM_AGE"

T0 <- ggplot(res_T0,
                  aes(x = logFC,
                      y = -log10(adj.P.Val),
                      label = delabel,
                      fill = diffexpressed,
                      shape = Direction_IMM_AGE)) +
  geom_point(size = 5, colour = "black") +
  scale_fill_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue",
    "NO" = "grey80"
  )) +
  scale_shape_manual(values = c(
    "> IMM_AGE" = 24,    
    "< IMM_AGE" = 21,   
    "Other" = 21
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text_repel(
    size = 8,
    force = 2,
    max.overlaps = Inf,        # allow many labels
    box.padding = 1.2,
    point.padding = 0.9,
    min.segment.length = 0,    # ALWAYS draw line if moved
    segment.color = "black",
    segment.size = 0.5
  ) +
  theme_classic() +
  xlim(-2.5, 2.5) +
  ylim(0, 5) +
  ggtitle("IMM-AGE Associated Genes \nPCV13 T0 v Placebo T0")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.position = "none"
    # legend.title = element_text(size = 20, face = "bold"),
    # legend.text  = element_text(size = 20, face = "bold"),
    # legend.key.size = unit(1.2, "cm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 6))
  ) +
  theme(plot.title = element_text(hjust=0.5, face = "bold"))

T0

###########################  Placebo T1 vs Placebo T0 ##########################

res_Placebo <- topTable(fit2, coef = "T1_Placebo_vs_T0_Placebo", number = Inf) %>% rownames_to_column("Gene_name")
# add a column of NAs
res_Placebo$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_Placebo$diffexpressed[res_Placebo$logFC > 0.5 & res_Placebo$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_Placebo$diffexpressed[res_Placebo$logFC < -0.5 & res_Placebo$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res_Placebo$delabel <- NA
res_Placebo$delabel[res_Placebo$diffexpressed != "NO"] <- res_Placebo$Gene_name[res_Placebo$diffexpressed != "NO"]

res_Placebo$Direction_IMM_AGE <- "Other"
res_Placebo$Direction_IMM_AGE[res_Placebo$Gene_name %in% genes_pos_fixed] <- "> IMM_AGE"
res_Placebo$Direction_IMM_AGE[res_Placebo$Gene_name %in% genes_neg_fixed] <- "< IMM_AGE"

Placebo <- ggplot(res_Placebo,
                  aes(x = logFC,
                      y = -log10(adj.P.Val),
                      label = delabel,
                      fill = diffexpressed,
                      shape = Direction_IMM_AGE)) +
  geom_point(size = 5, colour = "black") +
  scale_fill_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue",
    "NO" = "grey80"
  )) +
  scale_shape_manual(values = c(
    "> IMM_AGE" = 24,    
    "< IMM_AGE" = 21,   
    "Other" = 21
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text_repel(
    size = 8,
    force = 2,
    max.overlaps = Inf,        # allow many labels
    box.padding = 1.2,
    point.padding = 0.9,
    min.segment.length = 0,    # ALWAYS draw line if moved
    segment.color = "black",
    segment.size = 0.5
  ) +
  theme_classic() +
  xlim(-2.5, 2.5) +
  ylim(0, 5) +
  ggtitle("IMM-AGE Associated Genes \nPlacebo T1 v Placebo T0")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.position = "none"
    # legend.title = element_text(size = 20, face = "bold"),
    # legend.text  = element_text(size = 20, face = "bold"),
    # legend.key.size = unit(1.2, "cm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 6))
  ) +
  theme(plot.title = element_text(hjust=0.5, face = "bold"))

Placebo

###########################  PCV13 T1 vs PCV13 T0 ##############################

res_PCV13 <- topTable(fit2, coef = "T1_PCV13_vs_T0_PCV13", number = Inf) %>% rownames_to_column("Gene_name")
# add a column of NAs
res_PCV13$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_PCV13$diffexpressed[res_PCV13$logFC > 0.5 & res_PCV13$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_PCV13$diffexpressed[res_PCV13$logFC < -0.5 & res_PCV13$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res_PCV13$delabel <- NA
res_PCV13$delabel[res_PCV13$diffexpressed != "NO"] <- res_PCV13$Gene_name[res_PCV13$diffexpressed != "NO"]

res_PCV13$Direction_IMM_AGE <- "Other"
res_PCV13$Direction_IMM_AGE[res_PCV13$Gene_name %in% genes_pos_fixed] <- "> IMM_AGE"
res_PCV13$Direction_IMM_AGE[res_PCV13$Gene_name %in% genes_neg_fixed] <- "< IMM_AGE"

PCV13 <- ggplot(res_PCV13,
                aes(x = logFC,
                    y = -log10(adj.P.Val),
                    label = delabel,
                    fill = diffexpressed,
                    shape = Direction_IMM_AGE)) +
  geom_point(size = 5, colour = "black") +
  scale_fill_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue",
    "NO" = "grey80"
  )) +
  scale_shape_manual(values = c(
    "> IMM_AGE" = 24,     # circle
    "< IMM_AGE" = 21,     # triangle
    "Other" = 21
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text_repel(
    size = 8,
    force = 2,
    max.overlaps = Inf,        # allow many labels
    box.padding = 1.2,
    point.padding = 0.9,
    min.segment.length = 0,    # ALWAYS draw line if moved
    segment.color = "black",
    segment.size = 0.5
  ) +
  theme_classic() +
  xlim(-2.5, 2.5) +
  ylim(0, 5) +
  ggtitle("IMM-AGE Associated Genes \nPCV13 T1 v PCV13 T0")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.position = "none"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 6))
  ) +
  theme(plot.title = element_text(hjust=0.5, face = "bold"))

PCV13


################################################################################
#SAS and MAS 

#add in gene sets
SAS <- c("CCL4L2",
         "CCR4",
         "CCR7",
         "CD27",
         "CD40LG",
         "CXCL8",
         "CXCR5",
         "ETS1",
         "GPR183", 
         "HLA-DQA1",
         "HLA-DRB1",
         "HLA-DRB5",
         "ICOS",
         "IL24",
         "IL7R",
         "MS4A2",
         "PTGDR2",
         "SUSD2",
         "TCF7", "TNFRSF25", "VPREB3")

MAS <- c("ADAM17",
         "ADM",
         "ANG",
         "C5AR1",
         "CAMP", 
         "CD36", 
         "DEFA3",
         "DEFA4",
         "DEFB1",
         "HAVCR2",
         "H2BC4",#"HIST1H2BC",
         "H2BC5", #HIST1H2BD",
         "H2BC7",#"HIST1H2BF",
         "H2BC8", #"HIST1H2BG",
         "H2BC12", #HIST1H2BK",
         "H2BC21", #"HIST2H2BE",
         "HMGB2",
         "MYD88",
         "RNASE3",
         "TBK1",
         "TLR2",
         "TNFSF8")

mat_SAS_MAS <- read.csv("new/revisions/mat_SAS_MAS.csv", row.names = 1)
head(mat_SAS_MAS)

expr <- as.matrix(mat_SAS_MAS)

#load in metadata
metadata <- read.csv("D:/RH/VACIRISS/anon_metadata.csv") #%>% dplyr::select("anon_sample_ID", "Cohort_time")
metadata <- metadata %>% filter(Timepoint == "T0" | Timepoint == "T1")
metadata <- metadata %>% dplyr::select("anon_sample_ID", "Cohort_time")

# ensure sample order matches metadata
expr <- expr[, metadata$anon_sample_ID]
# samples in metadata but not in expression matrix
setdiff(metadata$anon_sample_ID, colnames(expr))

# samples in expression matrix but not in metadata
setdiff(colnames(expr), metadata$anon_sample_ID)

# factors
metadata$Cohort_time <- factor(metadata$Cohort_time)
design <- model.matrix(~ 0 + Cohort_time, data = metadata)
colnames(design) <- sub("^Cohort_time", "", colnames(design))
colnames(design)

colnames(design)
fit <- lmFit(expr, design)
contrasts <- makeContrasts(
  T0_PCV13_vs_Placebo = PCV13_T0 - Placebo_T0,
  T1_PCV13_vs_Placebo = PCV13_T1 - Placebo_T1,
  T1_Placebo_vs_T0_Placebo = Placebo_T1 - Placebo_T0,
  T1_PCV13_vs_T0_PCV13 = PCV13_T1 - PCV13_T0,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
res_T0 <- topTable(fit2, coef = "T0_PCV13_vs_Placebo", number = Inf)
res_T1 <- topTable(fit2, coef = "T1_PCV13_vs_Placebo", number = Inf)
res_Placebo <- topTable(fit2, coef = "T1_Placebo_vs_T0_Placebo", number = Inf)
res_PCV13 <- topTable(fit2, coef = "T1_PCV13_vs_T0_PCV13", number = Inf)

###########################  PCV13 T0 vs Placebo T0 ############################

res_T0 <- topTable(fit2, coef = "T0_PCV13_vs_Placebo", number = Inf) %>% rownames_to_column("Gene_name")
# add a column of NAs
res_T0$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_T0$diffexpressed[res_T0$logFC > 0.5 & res_T0$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_T0$diffexpressed[res_T0$logFC < -0.5 & res_T0$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res_T0$delabel <- NA
res_T0$delabel[res_T0$diffexpressed != "NO"] <- res_T0$Gene_name[res_T0$diffexpressed != "NO"]

res_T0$Signature <- "Other"
res_T0$Signature[res_T0$Gene_name %in% SAS] <- "SAS"
res_T0$Signature[res_T0$Gene_name %in% MAS] <- "MAS"

T0 <- ggplot(res_T0,
                  aes(x = logFC,
                      y = -log10(adj.P.Val),
                      label = delabel,
                      fill = diffexpressed,
                      shape = Signature)) +
  geom_point(size = 5, colour = "black") +
  scale_fill_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue",
    "NO" = "grey80"
  )) +
  scale_shape_manual(values = c(
    "SAS" = 21,     # circle
    "MAS" = 24,     # triangle
    "Other" = 21
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text_repel(
    size = 8,
    force = 2,
    max.overlaps = Inf,        # allow many labels
    box.padding = 1.2,
    point.padding = 0.9,
    min.segment.length = 0,    # ALWAYS draw line if moved
    segment.color = "black",
    segment.size = 0.5
  ) +
  
  theme_classic() +
  xlim(-2.5, 2.5) +
  ylim(0, 5) +
  ggtitle("SAS and MAS Associated Genes \nPCV13 T0 v Placebo T0")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.position = "none"
    # legend.title = element_text(size = 20, face = "bold"),
    # legend.text  = element_text(size = 20, face = "bold"),
    # legend.key.size = unit(1.2, "cm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 6))
  ) +
  theme(plot.title = element_text(hjust=0.5, face = "bold"))

T0

###########################  Placebo T1 vs Placebo T0 ##########################

res_Placebo <- topTable(fit2, coef = "T1_Placebo_vs_T0_Placebo", number = Inf) %>% rownames_to_column("Gene_name")
# add a column of NAs
res_Placebo$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_Placebo$diffexpressed[res_Placebo$logFC > 0.5 & res_Placebo$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_Placebo$diffexpressed[res_Placebo$logFC < -0.5 & res_Placebo$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res_Placebo$delabel <- NA
res_Placebo$delabel[res_Placebo$diffexpressed != "NO"] <- res_Placebo$Gene_name[res_Placebo$diffexpressed != "NO"]

res_Placebo$Signature <- "Other"
res_Placebo$Signature[res_Placebo$Gene_name %in% SAS] <- "SAS"
res_Placebo$Signature[res_Placebo$Gene_name %in% MAS] <- "MAS"

T0<- ggplot(res_Placebo,
             aes(x = logFC,
                 y = -log10(adj.P.Val),
                 label = delabel,
                 fill = diffexpressed,
                 shape = Signature)) +
  geom_point(size = 5, colour = "black") +
  scale_fill_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue",
    "NO" = "grey80"
  )) +
  scale_shape_manual(values = c(
    "SAS" = 21,     # circle
    "MAS" = 24,     # triangle
    "Other" = 21
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text_repel(
    size = 8,
    force = 2,
    max.overlaps = Inf,        # allow many labels
    box.padding = 1.2,
    point.padding = 0.9,
    min.segment.length = 0,    # ALWAYS draw line if moved
    segment.color = "black",
    segment.size = 0.5
  ) +
  
  theme_classic() +
  xlim(-2.5, 2.5) +
  ylim(0, 5) +
  ggtitle("SAS and MAS Associated Genes \nPlacebo T1 v Placebo T0")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.position = "none"
    # legend.title = element_text(size = 20, face = "bold"),
    # legend.text  = element_text(size = 20, face = "bold"),
    # legend.key.size = unit(1.2, "cm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 6))
  ) +
  theme(plot.title = element_text(hjust=0.5, face = "bold"))

T0

###########################  PCV13 T1 vs PCV13 T0 ##############################
res_PCV13 <- topTable(fit2, coef = "T1_PCV13_vs_T0_PCV13", number = Inf) %>% rownames_to_column("Gene_name")
# add a column of NAs
res_PCV13$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_PCV13$diffexpressed[res_PCV13$logFC > 0.5 & res_PCV13$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_PCV13$diffexpressed[res_PCV13$logFC < -0.5 & res_PCV13$adj.P.Val < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
res_PCV13$delabel <- NA
res_PCV13$delabel[res_PCV13$diffexpressed != "NO"] <- res_PCV13$Gene_name[res_PCV13$diffexpressed != "NO"]

res_PCV13$Signature <- "Other"
res_PCV13$Signature[res_PCV13$Gene_name %in% SAS] <- "SAS"
res_PCV13$Signature[res_PCV13$Gene_name %in% MAS] <- "MAS"

PCV13 <- ggplot(res_PCV13,
                aes(x = logFC,
                    y = -log10(adj.P.Val),
                    label = delabel,
                    fill = diffexpressed,
                    shape = Signature)) +
  
  geom_point(size = 5, colour = "black") +
  
  scale_fill_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue",
    "NO" = "grey80"
  )) +
  
  scale_shape_manual(values = c(
    "SAS" = 21,     # circle
    "MAS" = 24,     # triangle
    "Other" = 21
  )) +
  
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text_repel(
    size = 8,
    force = 2,
    max.overlaps = Inf,        # allow many labels
    box.padding = 1.2,
    point.padding = 0.9,
    min.segment.length = 0,    # ALWAYS draw line if moved
    segment.color = "black",
    segment.size = 0.5
  ) +
  
  theme_classic() +
  xlim(-2.5, 2.5) +
  ylim(0, 5) +
  ggtitle("SAS and MAS Associated Genes \nPCV13 T1 v PCV13 T0")+
  theme(
    plot.title   = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 18, face = "bold"),
    legend.position = "none"
    # legend.title = element_text(size = 20, face = "bold"),
    # legend.text  = element_text(size = 20, face = "bold"),
    # legend.key.size = unit(1.2, "cm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 6)),
    shape = guide_legend(override.aes = list(fill = "grey80", size = 6))
  ) +
  theme(plot.title = element_text(hjust=0.5, face = "bold"))

PCV13