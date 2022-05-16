#############################################
##                                         ##
##     DESeq and JTK Cycle on Males        ##
##                                         ##
#############################################

library(tidyverse)
library(MetaCycle)
library(DESeq2)
library(Bios2cor)
library(Biostrings)
library(tidysq)
library(data.table)
library(ggthemes)
library(plotly)

#############
#           #
#  DESeq2   #
#           #
#############

setwd("/Users/elizabethfreeman/Desktop/Mosquitoes/RNA-Seq")

males <- read_tsv("Male_Head_Cycling/male_head_counts.txt")

males_cleaned <- males

males_cleaned <- males_cleaned %>%
  mutate(males_cleaned, M_ZT0 = est_counts...2) %>%
  mutate(males_cleaned, M_ZT0.1 = est_counts...3) %>%
  mutate(males_cleaned, M_ZT0.2 = est_counts...4) %>%
  mutate(males_cleaned, M_ZT4 = est_counts...5) %>%
  mutate(males_cleaned, M_ZT4.1 = est_counts...6) %>%
  mutate(males_cleaned, M_ZT4.2 = est_counts...7) %>%
  mutate(males_cleaned, M_ZT8 = est_counts...8) %>%
  mutate(males_cleaned, M_ZT8.1 = est_counts...9) %>%
  mutate(males_cleaned, M_ZT8.2 = est_counts...10) %>%
  mutate(males_cleaned, M_ZT12 = est_counts...11) %>%
  mutate(males_cleaned, M_ZT12.1 = est_counts...12) %>%
  mutate(males_cleaned, M_ZT12.2 = est_counts...13) %>%
  mutate(males_cleaned, M_ZT16 = est_counts...14) %>%
  mutate(males_cleaned, M_ZT16.1 = est_counts...15) %>%
  mutate(males_cleaned, M_ZT16.2 = est_counts...16) %>%
  mutate(males_cleaned, M_ZT20 = est_counts...17) %>%
  mutate(males_cleaned, M_ZT20.1 = est_counts...18) %>%
  mutate(males_cleaned, M_ZT20.2 = est_counts...19)  %>%
  select(c(1, 20:37))

males_cleaned$M_ZT0 <- as.integer(males_cleaned$M_ZT0)
males_cleaned$M_ZT0.1 <- as.integer(males_cleaned$M_ZT0.1)
males_cleaned$M_ZT0.2 <- as.integer(males_cleaned$M_ZT0.2)
males_cleaned$M_ZT4 <- as.integer(males_cleaned$M_ZT4)
males_cleaned$M_ZT4.1 <- as.integer(males_cleaned$M_ZT4.1)
males_cleaned$M_ZT4.2 <- as.integer(males_cleaned$M_ZT4.2)
males_cleaned$M_ZT8 <- as.integer(males_cleaned$M_ZT8)
males_cleaned$M_ZT8.1 <- as.integer(males_cleaned$M_ZT8.1)
males_cleaned$M_ZT8.2 <- as.integer(males_cleaned$M_ZT8.2)
males_cleaned$M_ZT12 <- as.integer(males_cleaned$M_ZT12)
males_cleaned$M_ZT12.1 <- as.integer(males_cleaned$M_ZT12.1)
males_cleaned$M_ZT12.2 <- as.integer(males_cleaned$M_ZT12.2)
males_cleaned$M_ZT16 <- as.integer(males_cleaned$M_ZT16)
males_cleaned$M_ZT16.1 <- as.integer(males_cleaned$M_ZT16.1)
males_cleaned$M_ZT16.2 <- as.integer(males_cleaned$M_ZT16.2)
males_cleaned$M_ZT20 <- as.integer(males_cleaned$M_ZT20)
males_cleaned$M_ZT20.1 <- as.integer(males_cleaned$M_ZT20.1)
males_cleaned$M_ZT20.2 <- as.integer(males_cleaned$M_ZT20.2)

females <- read_tsv("Female_Head_Cycling/female_head_counts.txt")

females_cleaned <- females

females_cleaned <- females_cleaned %>%
  mutate(females_cleaned, F_ZT0 = est_counts...2) %>%
  mutate(females_cleaned, F_ZT0.1 = est_counts...3) %>%
  mutate(females_cleaned, F_ZT0.2 = est_counts...4) %>%
  mutate(females_cleaned, F_ZT4 = est_counts...5) %>%
  mutate(females_cleaned, F_ZT4.1 = est_counts...6) %>%
  mutate(females_cleaned, F_ZT4.2 = est_counts...7) %>%
  mutate(females_cleaned, F_ZT8 = est_counts...8) %>%
  mutate(females_cleaned, F_ZT8.1 = est_counts...9) %>%
  mutate(females_cleaned, F_ZT8.2 = est_counts...10) %>%
  mutate(females_cleaned, F_ZT12 = est_counts...11) %>%
  mutate(females_cleaned, F_ZT12.1 = est_counts...12) %>%
  mutate(females_cleaned, F_ZT12.2 = est_counts...13) %>%
  mutate(females_cleaned, F_ZT16 = est_counts...14) %>%
  mutate(females_cleaned, F_ZT16.1 = est_counts...15) %>%
  mutate(females_cleaned, F_ZT16.2 = est_counts...16) %>%
  mutate(females_cleaned, F_ZT20 = est_counts...17) %>%
  mutate(females_cleaned, F_ZT20.1 = est_counts...18) %>%
  mutate(females_cleaned, F_ZT20.2 = est_counts...19)  %>%
  select(c(1, 20:37))

females_cleaned$F_ZT0 <- as.integer(females_cleaned$F_ZT0)
females_cleaned$F_ZT0.1 <- as.integer(females_cleaned$F_ZT0.1)
females_cleaned$F_ZT0.2 <- as.integer(females_cleaned$F_ZT0.2)
females_cleaned$F_ZT4 <- as.integer(females_cleaned$F_ZT4)
females_cleaned$F_ZT4.1 <- as.integer(females_cleaned$F_ZT4.1)
females_cleaned$F_ZT4.2 <- as.integer(females_cleaned$F_ZT4.2)
females_cleaned$F_ZT8 <- as.integer(females_cleaned$F_ZT8)
females_cleaned$F_ZT8.1 <- as.integer(females_cleaned$F_ZT8.1)
females_cleaned$F_ZT8.2 <- as.integer(females_cleaned$F_ZT8.2)
females_cleaned$F_ZT12 <- as.integer(females_cleaned$F_ZT12)
females_cleaned$F_ZT12.1 <- as.integer(females_cleaned$F_ZT12.1)
females_cleaned$F_ZT12.2 <- as.integer(females_cleaned$F_ZT12.2)
females_cleaned$F_ZT16 <- as.integer(females_cleaned$F_ZT16)
females_cleaned$F_ZT16.1 <- as.integer(females_cleaned$F_ZT16.1)
females_cleaned$F_ZT16.2 <- as.integer(females_cleaned$F_ZT16.2)
females_cleaned$F_ZT20 <- as.integer(females_cleaned$F_ZT20)
females_cleaned$F_ZT20.1 <- as.integer(females_cleaned$F_ZT20.1)
females_cleaned$F_ZT20.2 <- as.integer(females_cleaned$F_ZT20.2)


both <- full_join(males_cleaned, females_cleaned, by = "target_id")


setwd("/Users/elizabethfreeman/Desktop/Mosquitoes/RNA-Seq/MvF_DESeq")

# Set up the conditions based on the experimental setup.

sex = rep(c('M','F'), c(18,18))
zt_m = rep(c('ZT0','ZT4','ZT8','ZT12','ZT16','ZT20'), c(3,3,3,3,3,3)) 
zt_f = rep(c('ZT0','ZT4','ZT8','ZT12','ZT16','ZT20'), c(3,3,3,3,3,3))
zt = c(zt_m, zt_f)

# Read the data from the standard input.
countData <- both %>% column_to_rownames("target_id")
  
# Build the dataframe from the conditions
samples = names(countData)
colData = data.frame(samples=samples, sex=sex, zt = zt)

#coldata = read.table(coldata_file, header=TRUE, sep="\t", row.names=1 )

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~ zt + sex)
#keep <- rowSums(counts(dds)) >=10
#dds <- dds[keep,]

#Set the reference to be compared
#dds$condition = relevel(dds$condition,"cond1")

# Run deseq
dds = DESeq(dds)

# Plot pca
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup=c("sex")) + geom_label(aes(label = rownames(colData))) #data point 30 is an outlier

#take out outlier and redo analysis

both <- both %>% column_to_rownames("target_id") %>% select(1:29, 31:36)

sex = rep(c('M','F'), c(18,17))
zt_m = rep(c('ZT0','ZT4','ZT8','ZT12','ZT16','ZT20'), c(3,3,3,3,3,3)) 
zt_f = rep(c('ZT0','ZT4','ZT8','ZT12','ZT16','ZT20'), c(3,3,3,2,3,3))
zt = c(zt_m, zt_f)

# Read the data from the standard input.
countData <- both 

# Build the dataframe from the conditions
samples = names(countData)
colData = data.frame(samples=samples, sex=sex, zt = zt)

#coldata = read.table(coldata_file, header=TRUE, sep="\t", row.names=1 )

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~ zt + sex)
#keep <- rowSums(counts(dds)) >=10
#dds <- dds[keep,]

#Set the reference to be compared
#dds$condition = relevel(dds$condition,"cond1")

# Run deseq
dds = DESeq(dds)

# Plot pca
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup=c("sex")) + geom_label(aes(label = rownames(colData))) #data point 12 is an outlier

# Format the results.
res = results(dds)

# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]

# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)
res05 = subset(sorted.df, sorted.df$padj <= 0.05)
res05.id = res05$id

# Write the table out.
write.table(sorted.df, file="ZTALL_BOTH_Head.new.deseq2.results.countsNotFiltered.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save the normalize data matrix.
write.table(dt, file="ZTALL_BOTH_Head.new.deseq2.normMatrix.countsNotFiltered.tsv", sep="\t",  row.names=FALSE,quote=FALSE)

# Select significant genes (res05.id genes) from normalized counts 
dt_sig = dt[dt$id %in% res05.id, ]
write.table(dt_sig, file="ZTALL_BOTH_Head.new.deseq2.normMatrix.sig.countsNotFiltered.tsv", sep="\t",  row.names=FALSE,quote=FALSE) 


#############
#           #
# Plotting  #
#           #
#############

setwd("/Users/elizabethfreeman/Desktop/Mosquitoes/RNA-Seq/")

vectorbase<- read_fasta("VectorBase-53_AgambiaePEST_AnnotatedTranscripts.fasta") #read in and clean genome
vb <- vectorbase %>%
  select(name) %>%
  separate(name, c("product", "gene", "organism", "gene_product", "transcript_product", 
                   "location", "length", "sequence", "SO", "is_pseudo"),  "\\|", extra = "merge")
vb$gene <- gsub(".*\\=", "", vb$gene)
vb$organism <- gsub(".*\\=", "", vb$organism)
vb$transcript_product <- gsub(".*\\=", "", vb$transcript_product)
vb$gene_product <- gsub(".*\\=", "", vb$gene_product)
vb$location <- gsub(".*\\=", "", vb$location)
vb$length <- gsub(".*\\=", "", vb$length)
vb$sequence <- gsub(".*\\=", "", vb$sequence)
vb$SO <- gsub(".*\\=", "", vb$SO)
vb$is_pseudo <- gsub(".*\\=", "", vb$is_pseudo)

vb <- vb %>% dplyr::rename(CycID = product)

vb$CycID <- str_trim(vb$CycID, side = "right")
vb$gene <- str_trim(vb$gene, side = "right")

goterms <- read_tsv("An_gambiae_GOTerms_250422.txt", col_names = FALSE) #read in GO Terms
goterms <- goterms %>% dplyr::rename(gene = X2, go = X5)
goterms <- goterms %>% select(c(gene, go)) #select for relevant columns

gochild <- read_csv("GO_Child_Terms.csv", col_names = FALSE)

gochild <- gochild %>% rename(go = X1)
match <- semi_join(goterms, gochild)
match <- match %>% group_by(gene) %>% nest()
match[,2] <- sapply(do.call(`c`,match[,2]), paste, collapse=", ")
match$data <- gsub("list", "", match$data)
match$data <- gsub("go = c\\(", "", match$data)
match$data <- gsub("\"", "", match$data)
match$data <- gsub("\\(", "", match$data)
match$data <- gsub("\\)", "", match$data)
match$data <- gsub("c", "", match$data)

match<- semi_join(vb, match)

#DESEq genes of interest
int_genes <- semi_join(dt, match, by = c("id" ="CycID"))

int_genes_sig <- semi_join(dt_sig, match, by = c("id" ="CycID"))

#combine and average timepoints 

ZTO <- tibble(int_genes_sig$F_ZT0, int_genes_sig$F_ZT0.1, int_genes_sig$F_ZT0.2)
int_genes_sig$ZT0_avg = base::rowMeans(ZTO)
ZT4 <- tibble(int_genes_sig$F_ZT4, int_genes_sig$F_ZT4.1, int_genes_sig$F_ZT4.2)
int_genes_sig$ZT4_avg = base::rowMeans(ZT4)
ZT8 <- tibble(int_genes_sig$F_ZT8, int_genes_sig$F_ZT8.1, int_genes_sig$F_ZT8.2)
int_genes_sig$ZT8_avg = base::rowMeans(ZT8)
ZT12 <- tibble(int_genes_sig$F_ZT12, int_genes_sig$F_ZT12.1)
int_genes_sig$ZT12_avg = base::rowMeans(ZT12)
ZT16 <- tibble(int_genes_sig$F_ZT16, int_genes_sig$F_ZT16.1, int_genes_sig$F_ZT16.2)
int_genes_sig$ZT16_avg = base::rowMeans(ZT16)
ZT20 <- tibble(int_genes_sig$F_ZT20, int_genes_sig$F_ZT20.1, int_genes_sig$F_ZT20.2)
int_genes_sig$ZT20_avg = base::rowMeans(ZT20)
ZTall <- tibble(int_genes_sig$ZT0_avg, int_genes_sig$ZT4_avg, int_genes_sig$ZT8_avg, 
                int_genes_sig$ZT12_avg, int_genes_sig$ZT16_avg, int_genes_sig$ZT20_avg)
int_genes_sig$ZTall_avg = base::rowMeans(ZTall)



long <- pivot_longer(int_genes_sig, cols = c(ZT0_avg, ZT4_avg, ZT8_avg, ZT12_avg, ZT16_avg, ZT20_avg),
                     names_to = "time", values_to = "avg")

ordered <- c("firebrick", "orange1", "darkgoldenrod3", "gold2", "forestgreen",
             "darkolivegreen4", "dodgerblue4", "deepskyblue1", 
             "cornflowerblue", "lightslateblue", "mediumpurple4", "lightpink", "thistle")

cycle_plot <- ggplot(data = long) +
  geom_line(aes(group = CycID, x = time, y = avg, color = CycID)) +
  scale_color_manual(values = ordered) +
  #scale_y_continuous(trans = "log10") +
  labs(x = "Z Time", y = "Log Averages of Normalized Counts") +
  theme_few()

ggplotly(cycle_plot)




write_csv(int_genes, "ZTALL_BOTH_Head_neurogenes.csv")

write_csv(int_genes_sig, "ZTALL_BOTH_Head_neurogenes_signficant.csv")








