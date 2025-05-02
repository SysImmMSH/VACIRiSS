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

#load metadata
metadata <- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata.csv") #400
metadata
colnames(metadata)


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
#txdb <- makeTxDbFromGFF("Transcript/gencode.v43.basic.annotation.gtf.gz", 
#                        format = "gtf")
#saveDb(txdb, file="Transcript/gencode.v43.sqlite")
#can just load this next time
txdb <- loadDb("Transcript/gencode.v43.sqlite")
#select columns required
columns(txdb)
k <- keys(txdb, "GENEID")
tx2gene <- select(txdb, keys = k, keytype = 'GENEID', columns = 'TXNAME')
head(tx2gene)

res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID")
tx2gene <- res[,2:1]
head(tx2gene)


#specify the working directory under the variable dir
dir <- "Transcript/transcript_level_quantification_anon/"
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

list.files("Transcript/transcript_level_quantification_anon/")

#use tximport to read our count files. Since we are using the same gtf file, the versions of the transcripts will be the same, hence the argument ignoreTxVersion = FALSE. If not, we set it to TRUE
txi_s <- tximport(files = files, type = 'salmon', tx2gene = tx2gene, importer=read.delim,
                  ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)

names(txi_s)
head(txi_s$counts)


# we set the column names to the sample names using our samples object

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

nrow(ddsTxi) #59735 genes measured

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
# Extract normalized counts
dds
dds1 <- dds

dds1

# save the normalised count data 
dds_df <- as.data.frame(counts(dds1,normalized = TRUE)) #%>% rownames_to_column("GENEID")
dds_df
#head(dds_df)
#dds_df$gene_ID <- substr(dds_df$GENEID,1,15) #need gene_id to allign with gene name

dds_df$gene_name <- mapIds(org.Hs.eg.db, rownames(dds_df), "SYMBOL", 
                           "ENSEMBL")
#dds_df<- dds_df %>% rownames_to_column("GENEID")

class(dds_df)
data<- dds_df
gene_ID <- rownames(dds_df)
rownames(data) <- NULL
data <- cbind(gene_ID,data)
data
data <- data %>% relocate(gene_name)

ncol(data)

write_csv(data, "paper/final anonymised dataframes/Transcript/Normalisedcounts_ddsDeSEQ2_Cohort_time_anon.csv")

########### some visualisation ##########
#this only works when the design only has one variable
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, c("Batch"))
plotPCA(vsd, c("Timepoint"))
plotPCA(vsd, c("Treatment"))
plotPCA(vsd, c("Cohort_time"))



library(limma)
library(edgeR)
#account for batch effect using limma (Matt did this on vsd object which used dds post DeSeq2)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Batch)
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
ggsave(plot = PCA_raw, "Figures/Transcript/New_transcript/PCA_limma_Cohort_time.png", width = 6, height = 5, dpi = 300)


#################### Start looking at differential expression, volcano  ##########
####using the dds object ####

resultsNames(dds)

levels(dds$Cohort_time)

##### control should be second option 
res <- results(dds, contrast = c("Cohort_time","PCV13_T0","Placebo_T0"))

res <- results(dds, contrast = c("Cohort_time","PCV13_T1","Placebo_T1"))

res <- results(dds, contrast = c("Cohort_time","PCV13_T2","Placebo_T2"))

res
res.sort <- res[order(res$pvalue),]
res.sort <- res[order(res$padj),]

res.sort 

summary(res)
head(res)
res
plotMA(res, ylim=c(-10,10))

plotCounts(dds, gene=which.min(res$padj), intgroup="Cohort_time")



write.csv(res, "Transcript/new_analysis/results_PCVT1vPlaceboT1_anon.csv")


#load annotation file
annotation <- read.csv("Transcript/annotation.csv")

timepoint <- read.csv("Transcript/new_analysis/results_PCVT1vPlaceboT1_anon.csv")

timepoint
#rename 
colnames(timepoint)
timepoint
names(timepoint)[names(timepoint) == "X"] <- "gene_ID"
colnames(annotation)

#join annotation file to ensembl id for significant DEGs
DEG_annot <- inner_join(annotation, timepoint, by = 'gene_ID')
summary(DEG_annot)
sum(is.na(DEG_annot$gene_name))
#replace blanks in gene_name with NA
DEG_annot[DEG_annot == ""] <- NA  
sum(is.na(DEG_annot$gene_name))

#get ensemble ID for NA genes
DEG_NA <- DEG_annot %>% 
  dplyr::filter(is.na (gene_name)) %>%
  dplyr::select(-gene_name) %>% 
  dplyr::rename(gene_name = 'gene_ID')

#join ensemble ID NAs with gene names and remove NAs from gene name
#DEG_annot <- full_join(DEG_NA, DEG_annot) %>% 
#  dplyr::filter(!is.na(gene_name))
#head(DEG_annot)

#removing 8th column (which is empty gene list)
#DEG_annot <- DEG_annot[, -8]
## try this way to add ensemble names to NA in gene name 


DEG_annot$gene_name <- ifelse(is.na(DEG_annot$gene_name), DEG_annot$gene_ID, DEG_annot$gene_name)

DEG_annot

write.csv(DEG_annot, "Transcript/new_analysis/results_PCVT1vPlaceboT1_anon.csv")

# select genes that are statistically significant FDR<0.05
DEG_p <- DEG_annot %>% 
  dplyr::filter(padj <0.05) %>%
  dplyr::select(gene_name, log2FoldChange) #1641 genes
summary(DEG_p)




######       Volcano Plot       ########
#####
####
###
##
#

#load df
df <- read.csv("Transcript/new_analysis/results_PCVT2vPlaceboT2_anon.csv")
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

VP_df <- ggplot(data=df, aes(x=(log2FoldChange), y=-log10(padj), label=delabel))+ 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 2)+
  theme_classic()+
  geom_vline(xintercept=c(-1, 1), linetype="dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", col="black")+
  scale_fill_manual(values=mycolors)+
  geom_label_repel(size = 2, force = 2, max.overlaps = 10)+
  ggtitle("PCV13 T2 v Placebo T2")+
  theme(plot.title = element_text(hjust=0.5, face = "bold"))+
  ylim(0,20)+
  xlim(-30, 30)

VP_df
ggsave(plot = VP_df, "Figures/Transcript/New_transcript/VP_PCVT0vPlaceboT0_limits_anon_new.png", dpi = 300, height = 10, width = 10, units = "in")

df
table(df$diffexpressed)

## sanity check for which way around 
d<-plotCounts(dds, gene="ENSG00000211955", intgroup="Cohort_time", returnData = TRUE)

d 

cyto_plot<-ggplot(d ,aes(x= Cohort_time, y= count)) +
  geom_boxplot(aes(fill= Cohort_time), width= 0.35, outlier.shape = NA, shape =21, alpha = 0.4)+
  geom_jitter(aes(fill= Cohort_time,  stroke = 1.5), width=0.25, shape = 21,alpha=1, size = 4) +
  #scale_color_manual(values= c("No"= "red", "Yes"= "green", "NA" = "grey" ))+
  scale_fill_manual(values= c("Placebo_T0" = "#bbcfc7",
                              "Placebo_T1" = "#85a89a",
                              "Placebo_T2" ="#577a6c",
                              "Placebo_T3" ="#30443c",
                              "PCV13_T0" = "#f3c6d9", 
                              "PCV13_T1" = "#e379a7",
                              "PCV13_T2" = "#d32c75",
                              "PCV13_T3" = "#ac2460"))+
  theme_classic()+
  #scale_y_log10()+ # limits= c(0.1, 10)
  #annotation_logticks(sides = "l")+
  ylab("Counts")+
  xlab("")+
  ggtitle("Test")+
  theme(plot.title = element_text(hjust=0.5))#+
#geom_vline(xintercept = which(WHO_T0_T2$Status == "Placebo T2"), linetype= "dotted", color = "black")
#geom_vline(xintercept = 4.5, linetype= "dotted", color = "black") +
#stat_compare_means(mapping= (aes(label = ..p.sig..)), comparisons = my_comps_test, method= "wilcox.test",  label = "p.signif", hide.ns = TRUE)
#geom_signif(comparisons = my_comps2, map_signif_level = TRUE, test = "wilcox.test")


cyto_plot
colnames(d)
# Calculate the average count per group
df_summary <- d %>%
  filter(!is.na(count)) %>%
  group_by(Cohort_time) %>%
  summarise(average_count = mean(count), .groups = "drop")
# View the result
df_summary


################# repeat for responders #####################


#load metadata
metadata <- read.csv("paper/final anonymised dataframes/Transcript/anon_metadata.csv") #400
metadata
table(metadata$Responder_final)

metadata_responders<- metadata %>%
  filter(!is.na(Responder_final))

# Make txdb object from GFF file
#txdb <- makeTxDbFromGFF("Transcript/gencode.v43.basic.annotation.gtf.gz", 
#                        format = "gtf")
#saveDb(txdb, file="Transcript/gencode.v43.sqlite")
#can just load this next time
txdb <- loadDb("Transcript/gencode.v43.sqlite")
#select columns required
columns(txdb)
k <- keys(txdb, "GENEID")
tx2gene <- select(txdb, keys = k, keytype = 'GENEID', columns = 'TXNAME')
head(tx2gene)

res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID")
tx2gene <- res[,2:1]
head(tx2gene)



#specify the working directory under the variable dir
dir <- "Transcript/transcript_level_quantification_anon/"
#create files (edit the metadata to different analyses - eg, day 0 only)
list.files(file.path(dir))
files <- file.path(dir, metadata_responders$name)
files
names(files) <- metadata_responders$anon_sample_ID
nrow(metadata_responders)
files

all(file.exists(files)) #should be true

files

list.files("Transcript/transcript_level_quantification_anon/")

#use tximport to read our count files. Since we are using the same gtf file, the versions of the transcripts will be the same, hence the argument ignoreTxVersion = FALSE. If not, we set it to TRUE
txi_s <- tximport(files = files, type = 'salmon', tx2gene = tx2gene, importer=read.delim,
                  ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)


names(txi_s)
head(txi_s$counts)

dim(metadata_responders)
# we now import the txi object into DESeq (~ tells DESEQ what to compare)


colnames(txi_s$counts) <- metadata_responders$anon_sample_ID

head(txi_s$counts)

metadata_responders$Cohort_time_resp<- paste(metadata_responders$Cohort_time, metadata_responders$Responder_final, sep = "_")
metadata_responders$Cohort_time_resp<- as.factor(metadata_responders$Cohort_time_resp)
metadata_responders$Cohort_time_resp

#### 
ddsTxi <- DESeqDataSetFromTximport(txi_s,
                                   colData = metadata_responders,
                                   design = ~ Batch  + Cohort_time_resp)



nrow(ddsTxi) #59735 genes measured


metadata_responders

#prefiltering (10 counts in at least 5 samples)
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 5]
nrow(ddsTxi)#23955 left after filtering


#run deseqobject
dds <- DESeq(ddsTxi)
resultsNames(dds)

levels(dds$Cohort_time)
levels(dds$Cohort_time_resp)

##############
dds1 <- dds


#remove final number from ensembl name so it matches with annotation file 
head(rownames(dds))
table(duplicated(substr(rownames(dds),1,15)))
rownames(dds) <- make.unique(substr(rownames(dds),1,15))
head(rownames(dds1))



#remove final number from ensembl name so it matches with annotation file 
head(rownames(dds1))
table(duplicated(substr(rownames(dds1),1,15)))
rownames(dds1) <- make.unique(substr(rownames(dds1),1,15))
head(rownames(dds1))

# save the normalised count data 
dds_df <- as.data.frame(counts(dds1,normalized = TRUE)) #%>% rownames_to_column("GENEID")
dds_df
#head(dds_df)
dds_df$gene_ID <- substr(dds_df$GENEID,1,15) #need gene_id to allign with gene name

dds_df$gene_name <- mapIds(org.Hs.eg.db, rownames(dds_df), "SYMBOL", 
                           "ENSEMBL")
#dds_df<- dds_df %>% rownames_to_column("GENEID")

class(dds_df)
data<- dds_df
gene_ID <- rownames(dds_df)
rownames(data) <- NULL
data <- cbind(gene_ID,data)
data
data <- data %>% relocate(gene_name)

write_csv(data, "paper/final anonymised dataframes/Transcript/Normalisedcounts_ddsDeSEQ2_responders_cohort_time_anon.csv")

####results(dds, contrast=c("condition","C","B")) meaning genes with logFC > 0 are overexpressed in C.  logFC < 0 means underexpressed in C and overexpressed in B

resultsNames(dds)

#res <- results(dds, name= "Cohort_time_resp_PCV13_T0_Yes_vs_PCV13_T0_No")
res <- results(dds,
               contrast = c("Cohort_time_resp", "PCV13_T0_Yes", "PCV13_T0_No"),
               lfcThreshold = 0, alpha = 0.1)

res <- results(dds,
               contrast = c("Cohort_time_resp", "PCV13_T1_Yes", "PCV13_T1_No"),
               lfcThreshold = 0, alpha = 0.1)

res <- results(dds,
               contrast = c("Cohort_time_resp", "PCV13_T2_Yes", "PCV13_T2_No"),
               lfcThreshold = 0, alpha = 0.1)


res
res.sort <- res[order(res$pvalue),]
res.sort <- res[order(res$padj),]


res.sort 

summary(res)
head(res)
res
plotMA(res, ylim=c(-10,10))





plotCounts(dds, gene="ENSG00000211813", intgroup="Cohort_time_resp")




####### change the gene names 


write.csv(res, "Transcript/new_analysis/results_Responders_cohort_time_T2_anon.csv")


#load annotation file
annotation <- read.csv("Transcript/annotation.csv")
timepoint <- read.csv("Transcript/new_analysis/results_Responders_cohort_time_T2_anon.csv")

timepoint
#rename 
colnames(timepoint)
timepoint
names(timepoint)[names(timepoint) == "X"] <- "gene_ID"
colnames(annotation)

#join annotation file to ensembl id for significant DEGs
DEG_annot <- inner_join(annotation, timepoint, by = 'gene_ID')
summary(DEG_annot)
sum(is.na(DEG_annot$gene_name))
#replace blanks in gene_name with NA
DEG_annot[DEG_annot == ""] <- NA  
sum(is.na(DEG_annot$gene_name))

#get ensemble ID for NA genes
DEG_NA <- DEG_annot %>% 
  dplyr::filter(is.na (gene_name)) %>%
  dplyr::select(-gene_name) %>% 
  dplyr::rename(gene_name = 'gene_ID')

#join ensemble ID NAs with gene names and remove NAs from gene name
#DEG_annot <- full_join(DEG_NA, DEG_annot) %>% 
#  dplyr::filter(!is.na(gene_name))
#head(DEG_annot)

#removing 8th column (which is empty gene list)
#DEG_annot <- DEG_annot[, -8]
## try this way to add ensemble names to NA in gene name 


DEG_annot$gene_name <- ifelse(is.na(DEG_annot$gene_name), DEG_annot$gene_ID, DEG_annot$gene_name)

DEG_annot

write.csv(DEG_annot, "Transcript/new_analysis/results_Responders_cohort_time_T2_anon.csv")

# select genes that are statistically significant FDR<0.05
DEG_p <- DEG_annot %>% 
  dplyr::filter(padj <0.05) %>%
  dplyr::select(gene_name, log2FoldChange) #1641 genes
summary(DEG_p)



######       Volcano Plot       ########
#####
####
###
##
#

#load df
df <- read.csv("Transcript/new_analysis/results_Responders_cohort_time_T1_anon.csv")
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

VP_df <- ggplot(data=df, aes(x=(log2FoldChange), y=-log10(padj), label=delabel))+ 
  geom_point(aes(fill = diffexpressed), pch = 21, colour = "black", size = 2)+
  theme_classic()+
  geom_vline(xintercept=c(-1, 1), linetype="dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", col="black")+
  scale_fill_manual(values=mycolors)+
  geom_label_repel(size = 2, force = 2, max.overlaps = 40)+
  ggtitle("PCV13 T0 Responder v Non-Responder")+
  theme(plot.title = element_text(hjust=0.5, face = "bold"))+
  ylim(0,50)+
  xlim(-50, 50)

VP_df
ggsave(plot = VP_df, "Figures/Transcript/New_transcript/VP_responder_T0_anon.png", dpi = 300, height = 10, width = 10, units = "in")

table(df$diffexpressed)
