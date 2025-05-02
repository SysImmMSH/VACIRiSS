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
library(ggplot2)
library(ggpubr)
library(biomaRt)

#load metadata
metadata <- read.csv("D:/RH/VACIRISS/anon_metadata.csv") #400
metadata

nrow(metadata)
table(metadata$Timepoint)

table(metadata$anon_sample_ID)


colnames(metadata)
# Check for duplicates in df
duplicated_rows_df <- metadata %>%
  group_by(anon_sample_ID) %>%
  filter(n() > 1) %>%
  ungroup()
duplicated_rows_df

# Make txdb object from GFF file
# txdb <- makeTxDbFromGFF("gencode.v43.basic.annotation.gtf.gz",
#                        format = "gtf")
# saveDb(txdb, file="gencode/gencode.v43.sqlite")

#can just load this next time
txdb <- loadDb("gencode/gencode.v43.sqlite")
#select columns required
columns(txdb)
k <- keys(txdb, "GENEID")
tx2gene <- select(txdb, keys = k, keytype = 'GENEID', columns = 'TXNAME')
head(tx2gene)

res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID")
tx2gene <- res[,2:1]
head(tx2gene)

#specify the working directory under the variable dir
dir <- "transcript_level_quantification_anon/transcript_level_quantification_anon"
#create files (edit the metadata to different analyses - eg, day 0 only)
list.files(file.path(dir))
files <- file.path(dir, metadata$name)
files
names(files) <- metadata$anon_sample_ID
nrow(metadata)
files

all(file.exists(files)) #should be true
#file.exists(files)
files

#use tximport to read our count files. Since we are using the same gtf file, the versions of the transcripts will be the same, hence the argument ignoreTxVersion = FALSE. If not, we set it to TRUE
txi_s <- tximport(files = files, type = 'salmon', tx2gene = tx2gene, importer=read.delim,
                  ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)

names(txi_s)
head(txi_s$counts)


#we set the column names to the sample names using our samples object
txi_s$counts
colnames(txi_s$counts) <- metadata$anon_sample_ID
head(txi_s$counts)
txi_s



metadata$Batch<- as.factor(metadata$Batch)
metadata$Cohort_time

# we now import the txi object into DESeq (~ tells DESEQ what to compare)
#design considers the effect of treatment/time etc

ddsTxi <- DESeqDataSetFromTximport(txi_s,
                                   colData = metadata,
                                   design = ~ Batch + Cohort_time)

nrow(ddsTxi)

metadata

#prefiltering (10 counts in at least 5 samples)
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 5]
nrow(ddsTxi)#26167 left after filtering


ddsTxi$Cohort_time <- relevel(ddsTxi$Cohort_time, ref = "Placebo_T0")


#run deseqobject
dds <- DESeq(ddsTxi)
resultsNames(dds)


head(rownames(dds))
dim(dds)

#remove final number from ensembl name so it matches with annotation file 
head(rownames(dds))
table(duplicated(substr(rownames(dds),1,15)))
rownames(dds) <- make.unique(substr(rownames(dds),1,15))
head(rownames(dds))

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, c("Batch"))
plotPCA(vsd, c("Timepoint"))
plotPCA(vsd, c("Treatment"))
plotPCA(vsd, c("Cohort_time"))

library(limma)
library(edgeR)
#account for batch effect using limma
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Batch)
mat <- assay(vsd)
mat

#nicer PCA - need to remove batch effect though
pcaData_raw <-plotPCA(vsd, intgroup=c("Treatment", "Timepoint", "Batch", "Cohort_time"), returnData=TRUE)
percentVar_raw <- round(100 * attr(pcaData_raw, "percentVar"))
#label_raw <- pcaData_raw %>% filter(case_id == "CC18" |
#                            case_id == "CC24" |
#                         case_id == "CC01" |
#                        case_id == "CC11" |
#                        case_id == "CC31"|
#                     case_id == "SC78")


PCA_raw = ggplot(pcaData_raw, aes(x = PC1, y =PC2))+
  geom_density_2d(colour = "black", alpha = 0.2)+
  geom_point(show.legend = TRUE, aes(fill = Cohort_time, shape = Treatment), size =3) +
  xlab(paste0("PC1: ",percentVar_raw[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_raw[2],"% variance")) +
  coord_fixed()+
  theme_bw()+
  #scale_fill_manual(values = c("1" = "cadetblue3",
  #                           "2" ="darkolivegreen1",
  #                          "3" = "darkmagenta"))+
  #scale_fill_manual(values = c("T0" = "#55a8ca",
  #                            "T1" ="orchid3",
  #                           "T2" = "#ede151",
  #                          "T3" = "mediumpurple"))+
  scale_shape_manual(values = c(21,22))+
  guides(fill = guide_legend(override.aes=list(shape=21))) #to colour the legend


PCA_raw
#ggsave(plot = PCA_raw, "Figures/Transcript/New_transcript/PCA_limma_Cohort_time.png", width = 6, height = 5, dpi = 300)


################################################################################

#preparing mat for gitools
colnames(mat)

df <- mat %>% as.data.frame()
head(df)

df$gene_name <- mapIds(org.Hs.eg.db, rownames(df), "SYMBOL", 
                       "ENSEMBL")

#add metadata
metadata <- read.csv("D:/RH/VACIRISS/anon_metadata.csv")
head(metadata)

table(metadata$Timepoint)
table(metadata$Cohort_time)

#make list of sample ids for baseline samples
sample_ids <- metadata %>%
  dplyr::filter(Timepoint == "T0") %>%
  dplyr::pull(anon_sample_ID)  # Extract sample IDs as a vector

# View the list of sample_ids
sample_ids

#extract baseline (T0) samples
df_T0<- df %>%
  dplyr::select(gene_name, all_of(sample_ids)) %>%
  filter(!is.na(gene_name))

df_T0


duplicates <- df_T0 %>%
  filter(!is.na(gene_name)) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup()

# Count duplicates
sum(duplicated(df_T0$gene_name))

#when gene_name is duplicated keep the row with most variation
df_T0$variance <- apply(df_T0[ , -1], 1, var)

# Keep row with max variance per gene
df_T0 <- df_T0 %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_max(order_by = variance, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(-variance)


df_T0_median_centered <- df_T0 %>%
  column_to_rownames("gene_name") %>%
  apply(1, function(x) {
    centered <- x - median(x, na.rm = TRUE)  # Subtract row median
    centered / sd(centered, na.rm = TRUE)    # Scale by row SD
  }) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("gene_name")

# View result
head(df_T0_median_centered)

#write.csv(df, "df.csv")
# write.table(df_T0_median_centered, "new/expression_data_T0_vst_batchcorrected_mc_sd.cdm", 
#             sep = "\t", row.names = FALSE, quote = FALSE)

library(cluster)
library(factoextra,  lib.loc = "D:/Win-Library/4.4")
library(NbClust,  lib.loc = "D:/Win-Library/4.4")

heatmap <- read.table("new/zscore_data_T0_vst_batchcorrected_mc_sd", header = TRUE, sep = "\t")
heatmap <- na.omit(heatmap)
head(heatmap)
rownames(heatmap) <- heatmap$X
heatmap <- heatmap[, -1] #%>% as.matrix()
heatmap 


set.seed(123)

mat <- t(heatmap)
fviz_nbclust(mat,
             FUNcluster = hcut,
             method = "gap_stat", 
             k.max = 10,         
             nboot = 500)

NbClust(mat, distance = "euclidean", min.nc=2, max.nc=10, 
        method = "ward.D", index = "all")
#optimum number of clusters = 2

distance_matrix <- dist(t(heatmap))
hc <- hclust(distance_matrix, method = "complete") 
plot(hc, main = "Dendrogram of Hierarchical Clustering", sub = "", xlab = "", ylab = "Height")
k <- 2
cluster_cut <- cutree(hc, k)
gr_df <- as.data.frame(cluster_cut) %>% rownames_to_column("sample_ID")
head(gr_df)
gr_df <- gr_df %>%
  mutate(cluster_cut = recode(cluster_cut, `1` = "B", `2` = "A"))
head(gr_df)

#png("dendrogram.png", width = 15, height = 8, units = "in", res = 300)  # High-resolution PNG
hc <- hclust(distance_matrix, method = "complete") 
plot(hc, main = "Dendrogram of Hierarchical Clustering", sub = "", xlab = "", ylab = "Height",  cex = 0.5)
rect.hclust(hc, k = 2, border = "blue")  # Cut into 3 clusters
#dev.off()  # Save and close the file
#dev.new()

################################################################################
#inflammation annotations
inflamGS <- c("HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_COMPLEMENT",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

#load in z-score output from gitools slea analysis
heatmap <- read.table("new/zscore_data_T0_vst_batchcorrected_mc_sd", header = TRUE, sep = "\t")
heatmap <- na.omit(heatmap)
head(heatmap)
rownames(heatmap) <- heatmap$X
heatmap <- heatmap[, -1] #%>% as.matrix()
heatmap 

inflamDF <- data.frame(median_inflam = apply(heatmap[inflamGS, ], MARGIN = 2, FUN = median)) %>% rownames_to_column()

global_median <- median(inflamDF$median_inflam)

head(inflamDF)

#to split into 2 groups
ha <- inflamDF %>%
  mutate(inflam.gs = findInterval(median_inflam,
                                  vec = quantile(median_inflam, probs = c(0, 0.5, 1)),
                                  rightmost.closed = TRUE),  # 2 groups: 1 (low), 2 (high)
         inflam.gs = ifelse(inflam.gs == 1, "low", "high")) %>%
  dplyr::select(rowname, inflam.gs)

head(ha)



#add to annotations
metadata <- read.csv("D:/RH/VACIRISS/anon_metadata.csv")
head(metadata)

responder <- read.csv("full_responder_df_anon_meta.csv") %>% dplyr::select("anon_sample_ID", "Responder_final")
head(responder)

table(metadata$Timepoint)
table(metadata$Cohort_time)

df_ann <- metadata %>%
  dplyr::filter(Timepoint == "T0") 
df_annotations <- df_ann[, c("anon_sample_ID", "Treatment", "Batch")]
df_annotations <- merge(df_annotations, ha, by.x= "anon_sample_ID", by.y = "rowname")
df_annotations <- merge(df_annotations, gr_df, by.x= "anon_sample_ID", by.y = "sample_ID") ###
df_annotations <- merge(df_annotations, responder, by = "anon_sample_ID")

df_annotations <- df_annotations %>% column_to_rownames("anon_sample_ID")

df_annotations <- df_annotations %>%
  rename_with(~ "Responders", matches("Responder_final")) %>%
  rename_with(~ "Clusters", matches("cluster_cut"))

df_annotations$Responders <- as.factor(df_annotations$Responders)
df_annotations$Treatment <- as.factor(df_annotations$Treatment)
df_annotations$Clusters <- as.factor(df_annotations$Clusters)
df_annotations <- df_annotations[, c("Treatment", "inflam.gs", "Clusters")]


head(df_annotations)
max_value <- max(abs(heatmap), na.rm = TRUE)  # Ensures symmetric color scaling
breaks <- seq(-max_value, max_value, length.out = 101)  

tier_colours <- c("low" = "#FFEB3B", "high" = "red")
treatment_colours <- c("PCV13" = "#d32c75", "Placebo" = "#577a6c")
cluster_colours <- c("A" = "#0dde3a", "B" = "#0003C7")

library(pheatmap)

#add * to   "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
# Define the list of row names to modify
hallmarks_to_star <- c(
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

# Add asterisk to matching row names
rownames(heatmap) <- sapply(rownames(heatmap), function(name) {
  if (name %in% hallmarks_to_star) {
    paste0(name, "*")
  } else {
    name
  }
})


#png("new/heatmap_zscore_slea_2cluster_vst_long_2inflam_tiers.png", width=25, height=35, units = "in", res = 500)
pheatmap((heatmap),
         fontsize_col = 20,  # Increase the column label font size
         fontsize_row = 22,  # Increase the row label font size
         fontsize = 20,
         main = "Baseline Endotypes",
         annotation_col = df_annotations,
         show_colnames = FALSE,
         color = colorRampPalette(c("#163f76", "#4575B4", "#FFFFF2", "#ec763f", "#bf2a1b"))(100),
         breaks = breaks,
         cluster_rows = TRUE,
         annotation_colors = list(inflam.gs = tier_colours, Treatment = treatment_colours, Clusters = cluster_colours),
         cutree_cols = 2)
#dev.off()

################################################################################
###working out proportions
head(df_annotations)

write.csv(df_annotations, "new/df_annotations_2tiers.csv", row.names = TRUE)


#proportion endotypes
df_annotations %>%
  dplyr::count(Clusters) %>%
  mutate(percentage = (n / sum(n)) * 100)

df_annotations %>%
  dplyr::count(Clusters) %>%
  mutate(percentage = round((n / sum(n)) * 100, 2)) %>%
  write.csv("new/cluster_cut_percentages_2tiers.csv", row.names = FALSE)

df_annotations %>%
  dplyr::count(inflam.gs) %>%
  mutate(percentage = round((n / sum(n)) * 100, 2)) %>%
  write.csv("new/inflam_percentages_2tiers.csv", row.names = FALSE)

#proportion inflammation tiers
df_annotations %>%
  group_by(Clusters, inflam.gs) %>%
  summarise(total = n(), .groups = 'drop') %>%
  group_by(Clusters) %>%
  mutate(proportion = round((total / sum(total)) * 100, 2)) %>%
  ungroup() %>%
  dplyr::select(Clusters, inflam.gs, total, proportion) %>%
  write.csv("new/tier_proportions_by_cluster_2tiers.csv", row.names = FALSE)


#treatment allocation by clusters and tiers
df_annotations %>%
  group_by(Clusters, Treatment) %>%
  summarise(total = n(), .groups = 'drop') %>%
  group_by(Clusters) %>%
  mutate(proportion = round((total / sum(total)) * 100, 2)) %>%
  ungroup() %>%
  dplyr::select(Clusters, Treatment, total, proportion) %>%
  write.csv("new/treatment_proportions_by_clusters.csv", row.names = FALSE)
