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
library(viridis)

#############
#           #
#  DESeq2   #
#           #
#############

setwd("/Users/elizabethfreeman/Desktop/Mosquitoes/RNA-Seq/Male_Head_Cycling")


males <- read_tsv("male_head_counts.txt")

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

countdata<- males_cleaned %>% column_to_rownames("target_id")

samples = "3x3x3x3x3x3"
# Extract the experimental design from the command line.
design = unlist(strsplit(samples, 'x'))

cond1_num = as.integer(design[1])
cond2_num = as.integer(design[2])
cond3_num = as.integer(design[3])
cond4_num = as.integer(design[4])
cond5_num = as.integer(design[5])
cond6_num = as.integer(design[6])

cond_1 = rep("ZT0", cond1_num)
cond_2 = rep("ZT4", cond2_num)
cond_3 = rep("ZT8", cond3_num)
cond_4 = rep("ZT12", cond4_num)
cond_5 = rep("ZT16", cond5_num)
cond_6 = rep("ZT20", cond6_num)

condition = factor(c(cond_1, cond_2, cond_3, cond_4, cond_5, cond_6))

# Build the dataframe from the conditions
samples = names(countdata)
colData = data.frame(samples=samples, condition=condition)

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countdata, colData=colData, design = ~condition)

dds <- DESeq(dds, test="LRT", reduced=~1)
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup=c("condition"))

# Format the results.
res = results(dds)

# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj)), ]

# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)
res05 = subset(sorted.df, sorted.df$padj <= 0.05)
res05.id = res05$id

# Write the table out.
write.table(sorted.df, file="ZTALL_Male_Head.new.deseq2.results.countsNotFiltered.tsv", sep="\t", row.names = FALSE, quote=FALSE)

# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)

# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)

# Save the normalize data matrix.
write.table(dt, file="ZTALL_Male_Head.new.deseq2.normMatrix.countsNotFiltered.tsv", sep="\t",  row.names=FALSE,quote=FALSE)

# Select significant genes (res05.id genes) from normalized counts 
dt_sig = dt[dt$id %in% res05.id, ]
write.table(dt_sig, file="ZTALL_Male_Head.new.deseq2.normMatrix.sig.tsv", sep="\t",  row.names=FALSE,quote=FALSE) 

#############
#           #
#   JTK     #
#           #
#############

timepoints = c(0,0,0,4,4,4,8,8,8,12,12,12,16,16,16,20,20,20)
meta2d(infile = "ZTALL_Male_Head.new.deseq2.normMatrix.sig.tsv", filestyle = "txt", outdir = "M_Head_MetaCycle", timepoints = timepoints, minper = 8, maxper=24, cycMethod = 'JTK')

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

match<- inner_join(vb, match)

M_JTK<- read_tsv("Male_Head_Cycling/M_Head_MetaCycle/JTKresult_ZTALL_Male_Head.new.deseq2.normMatrix.sig.tsv")

sigM_JTK <- subset(M_JTK, M_JTK$BH.Q <= 0.05) #only want significant reads

vb_m_JTK <- left_join(sigM_JTK, vb) #combine reads and genome


#DESEq genes of interest
int_genes <- semi_join(dt, match, by = c("id" ="CycID"))

#Metacycle genes of interest
int_cycle <- semi_join(vb_m_JTK, match)

int_cycle_nc <- left_join(int_cycle, dt, by = c("CycID" = "id"))

#combine and average timepoints 

ZTO <- tibble(int_cycle_nc$M_ZT0, int_cycle_nc$M_ZT0.1, int_cycle_nc$M_ZT0.2)
int_cycle_nc$ZT0_avg = base::rowMeans(ZTO)
ZT4 <- tibble(int_cycle_nc$M_ZT4, int_cycle_nc$M_ZT4.1, int_cycle_nc$M_ZT4.2)
int_cycle_nc$ZT4_avg = base::rowMeans(ZT4)
ZT8 <- tibble(int_cycle_nc$M_ZT8, int_cycle_nc$M_ZT8.1, int_cycle_nc$M_ZT8.2)
int_cycle_nc$ZT8_avg = base::rowMeans(ZT8)
ZT12 <- tibble(int_cycle_nc$M_ZT12, int_cycle_nc$M_ZT12.1, int_cycle_nc$M_ZT12.2)
int_cycle_nc$ZT12_avg = base::rowMeans(ZT12)
ZT16 <- tibble(int_cycle_nc$M_ZT16, int_cycle_nc$M_ZT16.1, int_cycle_nc$M_ZT16.2)
int_cycle_nc$ZT16_avg = base::rowMeans(ZT16)
ZT20 <- tibble(int_cycle_nc$M_ZT20, int_cycle_nc$M_ZT20.1, int_cycle_nc$M_ZT20.2)
int_cycle_nc$ZT20_avg = base::rowMeans(ZT20)
ZTall <- tibble(int_cycle_nc$ZT0_avg, int_cycle_nc$ZT4_avg, int_cycle_nc$ZT8_avg, 
                int_cycle_nc$ZT12_avg, int_cycle_nc$ZT16_avg, int_cycle_nc$ZT20_avg)
int_cycle_nc$ZTall_avg = base::rowMeans(ZTall)


write_csv(int_cycle_nc, "Male_Cycling_genesofinterest.csv")



long <- pivot_longer(int_cycle_nc, cols = c(ZT0_avg, ZT4_avg, ZT8_avg, ZT12_avg, ZT16_avg, ZT20_avg),
                     names_to = "time", values_to = "avg") %>%
  mutate(time = factor(time, levels = c("ZT0_avg", "ZT4_avg", "ZT8_avg", "ZT12_avg",
                                        "ZT16_avg", "ZT20_avg")))


ordered <- c("firebrick", "indianred4", "indianred2", "coral3", "coral", "orange4", "orange1", "darkgoldenrod3", "gold2", "khaki4", "forestgreen",
             "darkolivegreen4", "darkolivegreen", "dodgerblue4", "dodgerblue", "deepskyblue4", "deepskyblue1", 
             "cornflowerblue", "lightslateblue", "mediumpurple4", "mediumpurple2","orchid4", "orchid", "lightpink", "thistle", "thistle2",
             "tan4", "burlywood4", "burlywood1")
cycle_plot <- ggplot(data = long) +
  geom_line(aes(group = CycID, x = time, y = avg, color = CycID)) +
  scale_color_manual(values = ordered) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Z Time", y = "Log Averages of Normalized Counts") +
  theme_few()

ggplotly(cycle_plot)

int_genes_sig <- semi_join(dt_sig, match, by = c("id" ="CycID"))

ZTO <- tibble(int_genes_sig$M_ZT0, int_genes_sig$M_ZT0.1, int_genes_sig$M_ZT0.2)
int_genes_sig$ZT0_avg = base::rowMeans(ZTO)
ZT4 <- tibble(int_genes_sig$M_ZT4, int_genes_sig$M_ZT4.1, int_genes_sig$M_ZT4.2)
int_genes_sig$ZT4_avg = base::rowMeans(ZT4)
ZT8 <- tibble(int_genes_sig$M_ZT8, int_genes_sig$M_ZT8.1, int_genes_sig$M_ZT8.2)
int_genes_sig$ZT8_avg = base::rowMeans(ZT8)
ZT12 <- tibble(int_genes_sig$M_ZT12, int_genes_sig$M_ZT12.1, int_genes_sig$M_ZT12.2)
int_genes_sig$ZT12_avg = base::rowMeans(ZT12)
ZT16 <- tibble(int_genes_sig$M_ZT16, int_genes_sig$M_ZT16.1, int_genes_sig$M_ZT16.2)
int_genes_sig$ZT16_avg = base::rowMeans(ZT16)
ZT20 <- tibble(int_genes_sig$M_ZT20, int_genes_sig$M_ZT20.1, int_genes_sig$M_ZT20.2)
int_genes_sig$ZT20_avg = base::rowMeans(ZT20)
ZTall <- tibble(int_genes_sig$ZT0_avg, int_genes_sig$ZT4_avg, int_genes_sig$ZT8_avg, 
                int_genes_sig$ZT12_avg, int_genes_sig$ZT16_avg, int_genes_sig$ZT20_avg)
int_genes_sig$ZTall_avg = base::rowMeans(ZTall)

long_genes <- pivot_longer(int_genes_sig, cols = c(ZT0_avg, ZT4_avg, ZT8_avg, ZT12_avg, ZT16_avg, ZT20_avg),
                           names_to = "time", values_to = "avg") %>%
  mutate(time = factor(time, levels = c("ZT0_avg", "ZT4_avg", "ZT8_avg", "ZT12_avg",
                                        "ZT16_avg", "ZT20_avg")))

gene_plot <- ggplot(data = long_genes) +
  geom_line(aes(group = id, x = time, y = avg, color = id)) +
 # scale_color_manual(values = ordered) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Z Time", y = "Log Averages of Normalized Counts") +
  theme_few()

ggplotly(gene_plot)


write_csv(int_genes, "ZTALL_Male_Head_neurogenes.csv")
write_csv(int_genes_sig, "ZTALL_Male_Head_neurogenes_signficant.csv")
write_csv(int_cycle_nc, "ZTALL_Male_Head_neurogenes_cycling.csv")
dt_info <- full_join(int_genes, match, by = c("id" ="CycID"))
write_csv(dt_info, "ZTALL_Male_Head_neurogenes.csv")




#Normalized expression

#this is the avg of the times/ maximum count within that time 

#let's do it with cycling 

int_cycle_nc$NEZT0 <- int_cycle_nc$ZT0_avg / pmax(int_cycle_nc$M_ZT0, int_cycle_nc$M_ZT0.1, int_cycle_nc$M_ZT0.2)
int_cycle_nc$NEZT4 <- int_cycle_nc$ZT4_avg / pmax(int_cycle_nc$M_ZT4, int_cycle_nc$M_ZT4.1, int_cycle_nc$M_ZT4.2)
int_cycle_nc$NEZT8 <- int_cycle_nc$ZT8_avg / pmax(int_cycle_nc$M_ZT8, int_cycle_nc$M_ZT8.1, int_cycle_nc$M_ZT8.2)
int_cycle_nc$NEZT12 <- int_cycle_nc$ZT12_avg / pmax(int_cycle_nc$M_ZT12, int_cycle_nc$M_ZT12.1, int_cycle_nc$M_ZT12.2)
int_cycle_nc$NEZT16 <- int_cycle_nc$ZT16_avg / pmax(int_cycle_nc$M_ZT16, int_cycle_nc$M_ZT16.1, int_cycle_nc$M_ZT16.2)
int_cycle_nc$NEZT20 <- int_cycle_nc$ZT20_avg / pmax(int_cycle_nc$M_ZT20, int_cycle_nc$M_ZT20.1, int_cycle_nc$M_ZT20.2)

long_NE_cycling <- pivot_longer(int_cycle_nc, cols = c(NEZT0 , NEZT4, NEZT8, NEZT12, NEZT16, NEZT20),
                                names_to = "time", values_to = "avg") %>%
  mutate(time = factor(time, levels = c("NEZT0", "NEZT4", "NEZT8", "NEZT12",
                                        "NEZT16", "NEZT20")))

ggplot(data = long_NE_cycling) + 
  geom_tile(aes(x = time, y = CycID, fill = avg)) +
  scale_fill_viridis(option = "plasma") +
  theme_few()

