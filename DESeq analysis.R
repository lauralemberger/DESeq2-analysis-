install.packages("dplyr")
library(dplyr)


#installing DESeq2 package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
force = TRUE
browseVignettes("DESeq2")
install.packages("DESeq2", type = "source")

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

install.packages("tidyverse")

library(DESeq2)

install.packages("devtools")
devtools::install_github("r-lib/conflicted")

library(tidyverse)

install.packages("airway")
library(airway)

#read in data

d <- read.delim("multibamsummary1.csv", sep = ";", header = TRUE)
counts_data <- paste0(d$X..chr, ":", d$X.start., "-", d$X.end.)
row.names(d) <- counts_data
d <- d %>% select(-X..chr, -X.start., -X.end.)
print(d)
dim(d)
colnames(d)
d
colData <- read.csv2("sample_info.csv", header = TRUE)

print(d)
print(colData)
rownames(colData) <- colData$X
colnames(d) == rownames(colData)
colnames(d)[1]
rownames(colData)[1]

str_split_i(colnames(d),"\\.",2)
colnames(d) <- str_split_i(colnames(d),"\\.",2)

str_split_i(rownames(colData),"\\.",2)
rownames(colData) <- str_split_i(rownames(colData),"\\.",2)



#construct the deseqdataset

dds <- DESeqDataSetFromMatrix(countData = d,
                              colData = colData,
                              design = ~ construct + replica + time)

dds
colData(dds)
head(assay(dds))
head(d)

#rlog transformation 

rlog(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

#PCA Plot

plotPCA(rld, intgroup = c("construct", "time"))
pcaData <- plotPCA(rld, intgroup = c( "construct", "time"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, main = "PCA", aes(x = PC1, y = PC2, color = as.factor(time), shape = construct)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rlog data")

#heatmap sample-to-sample-distances

install.packages("pheatmap")

library("pheatmap")

library("RColorBrewer")

sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$construct, rld$time, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#likelihood ratio test (Hypothesentest)

datasetdeseq <- DESeq(dds, test="LRT", reduced=~construct+replica)

res <- results(datasetdeseq)
res #beinhaltet alle Peaks 
summary(res, alpha=0.05)
view(res)
res[c(1,5),]

#filtering significant peaks

sum(res$padj < 0.05, na.rm=TRUE)

data.significant <- subset(res, padj < 0.05)
summary(data.significant)
print(data.significant)
head(data.significant)
print(data.significant)
data.significant

#visualization of significant peaks
head(assay(rld))
sortCol<-colnames(rld)[order(rld$construct, rld$time)]

df <- as.data.frame(colData(dds)[,c("construct","time")])
pheatmap(assay(rld)[rownames(data.significant),sortCol], show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, scale = "row")

length(select)


#Tabelle mit bekannten Genomregionen (signifikante Signale/Peaks) in IGV visualisieren
#dafür zuerst die Tabelle in ein bed file umwandeln 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(rtracklayer)

#rownames filtern in Vektoren: Chr, Start, End

print(data.significant)
rownames(data.significant) 
significantData <- rownames(data.significant)
view(significantData)

str_split(significantData, ":" )
str_split_i(significantData, ":", 1)
chr <- str_split_i(significantData, ":", 1)

str_split_i(significantData, ":", 2)
region <- str_split_i(significantData, ":", 2)
str_split(region, "-")
str_split_i(region, "-", 1)
start <- str_split_i(region, "-", 1)
view(start)

str_split(significantData, "-")
str_split_i(significantData, "-", 2)                                                      
end <- str_split_i(significantData, "-", 2)

x <- GRanges(chr, ranges = IRanges(start = as.numeric(start),
                                   end = as.numeric(end), 
                                   score = data.significant$padj))
x
export.bed(x, con = "C:/Users/laura/Documents/Bioinformatik Praktikum/significantpeaks.bed")
save(x, file =  "C:/Users/laura/Documents/Bioinformatik Praktikum/x.rda")



#plotCounts 

time_factor <- colnames(colData)
time_factor <- paste(datasetdeseq$time, datasetdeseq$construct)
time_factor
str_split_i(time_factor, " ", 1)
timex <- str_split_i(time_factor, " ", 1)
timex
timex <- factor(timex, labels = c("1", "2", "0", "4"))
timex

timex <- as.factor(timex)
typeof(datasetdeseq$time)
typeof(timex)
class(timex)

colData(datasetdeseq)

plotCounts(datasetdeseq, gene = "chr18:62420542-62421254", 
           intgroup = "time")

time_factor2 <- paste(datasetdeseq$construct, datasetdeseq$time, sep = "_")
time_factor2

datasetdeseq$time_construct <- as.factor(time_factor2)
colData(datasetdeseq)

#Regionen aus significantData nehmen!

plotCounts(datasetdeseq, gene = "chr14:64840548-64843743",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr11:8831934-8834702",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr20:48846149-48847290",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr1:116508159-116509891",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr22:33919379-33921043",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr14:34718145-34719716",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr2:157874505-157877209",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr3:29683117-29684597",
           intgroup = "time_construct")

plotCounts(datasetdeseq, gene = "chr6:118934486-118935565", 
           intgroup = "time_construct")

#rausfiltern der kleinsten p-Values, um eindeutigere 
#Beispiele zu bekommen (Histonacetylierung nimmt zu / ab)

smallest_p_values <- subset(res, padj < 0.01)
smallest_p_values
view(smallest_p_values)

#Zeitpunkt 0, dann Abnahme, dann wieder Zunahme bis 4
plotCounts(datasetdeseq, gene = "chr5:37667044-37670097", 
           intgroup = "time_construct")

#Abnahme Histon Acetylierung
plotCounts(datasetdeseq, gene = "chr20:43666580-43668593", 
           intgroup = "time_construct")

#Zunahme Histon Acetylierung
plotCounts(datasetdeseq, gene = "chr1:47407522-47412846", 
           intgroup = "time_construct")

#Abnahme
plotCounts(datasetdeseq, gene = "chr16:10331550-10337729", 
           intgroup = "time_construct")


smallest_p_values2 <- subset(res, padj < 0.0001)
smallest_p_values2
view(smallest_p_values2)

#Zunahme
plotCounts(datasetdeseq, gene = "chr14:68783769-68797542", 
           intgroup = "time_construct")

#Abnahme
plotCounts(datasetdeseq, gene = "chr5:88784399-88788039", 
           intgroup = "time_construct")

#Abnahme
plotCounts(datasetdeseq, gene = "chr19:50003227-50005416", 
           intgroup = "time_construct")

smallest_p_values3 <- subset(res, padj < 0.000001)
smallest_p_values3
view(smallest_p_values3)

#Abnahme
plotCounts(datasetdeseq, gene = "chr3:191328323-191332421", 
           intgroup = "time_construct")

#Zunahme
plotCounts(datasetdeseq, gene = "chr15:60399256-60403675", 
           intgroup = "time_construct")

#Abnahme
plotCounts(datasetdeseq, gene = "chr4:48339542-48343224", 
           intgroup = "time_construct")

#Zunahme
plotCounts(datasetdeseq, gene = "chr1:149843256-149850557", 
           intgroup = "time_construct")

#Zunahme
plotCounts(datasetdeseq, gene = "chr8:127734251-127741103", 
           intgroup = "time_construct")

