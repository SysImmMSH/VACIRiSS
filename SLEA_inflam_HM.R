#SLEA heatmap with inflammatory subtypes
library(cluster)
library(factoextra)
library(NbClust)


#clustering
heatmap <- read.csv("new/zscore_data_Gitools.csv")
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

################################################################################
#visualising SLEA results in heatmap
#inflammation annotations
inflamGS <- c("HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_COMPLEMENT",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING",
              "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

#load in z-score output from gitools slea analysis
heatmap <- read.csv("new/zscore_data_Gitools.csv")
heatmap <- na.omit(heatmap)
head(heatmap)
rownames(heatmap) <- heatmap$X
heatmap <- heatmap[, -1]
heatmap 

inflamDF <- data.frame(median_inflam = apply(heatmap[inflamGS, ], MARGIN = 2, FUN = median)) %>% rownames_to_column()

global_median <- median(inflamDF$median_inflam)

head(inflamDF)

#split into 2 groups
ha <- inflamDF %>%
  mutate(inflam.gs = findInterval(median_inflam,
                                  vec = quantile(median_inflam, probs = c(0, 0.5, 1)),
                                  rightmost.closed = TRUE),  # 2 groups: 1 (low), 2 (high)
         inflam.gs = ifelse(inflam.gs == 1, "low", "high")) %>%
  dplyr::select(rowname, inflam.gs)

head(ha)

#add to annotations
metadata <- read.csv("anon_metadata.csv")
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
#Define the list of row names to modify
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


p <- pheatmap((heatmap),
              fontsize_col = 18,
              fontsize_row = 20,
              fontsize = 18,
              main = "Baseline Endotypes",
              annotation_col = df_annotations,
              show_colnames = FALSE,
              color = colorRampPalette(c("#163f76", "#4575B4", "#FFFFF2", "#ec763f", "#bf2a1b"))(100),
              breaks = breaks,
              cluster_rows = TRUE,
              annotation_colors = list(inflam.gs = tier_colours, Treatment = treatment_colours, Clusters = cluster_colours),
              cutree_cols = 2)

