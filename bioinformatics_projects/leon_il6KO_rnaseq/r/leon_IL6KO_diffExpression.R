#load library in
library(edgeR)
library(ggplot2)

#load in metadata
#add short names as well as file names
samples= read.table("/Users/davehill/Desktop/LEON_01222022_RNASEQ_IL6KO/counts_updated/LEON_METADATA.txt", sep="\t", header=TRUE)

#creates a vector of sample names
sampleNames <- samples$Sample

#load in gene length data
lengths <- read.csv("/Users/davehill/Desktop/gene_lengths/GRCm39_lengths.txt", sep="\t", header=TRUE)

#read count files into DGE object
counts = readDGE(samples$countFile)

# indicies of metatags
noint = rownames(counts) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")

# counts per million to exclude low expressers
cpms = cpm(counts)

# keep rows that are not metadata tags and have more expressions than 1
keep = rowSums(cpms >1) >=3 & !noint

#filtered counts matrix
counts = counts[keep, ]

#make DGElist to hold expression data
d <- DGEList(counts=counts, group = samples$Treatment)

#Adds the short names to the d object
d$samples$name <- samples$ShortName

#add data frame to include gene names and their lengths( for RPKM function)
d$genes <- data.frame(gene=rownames(d), length=lengths$length[match(rownames(d), lengths$gene)])

#estimate normalization factors
d <- calcNormFactors(d)

#MDS plot
plotMDS(d,pch=1, labels=samples$ShortName,main = "MDS for top 500 genes", col=c("blue", "red", "purple", "green", "orange", "yellow", "violet")[factor(samples$Treatment)], xlim=range(-2.5:2.5), ylim=range(-2.5:2.5), top=500)
#plotMDS(d, labels=samples$ShortName, col=c("blue", "red", "purple", "green", "orange", "yellow", "violet")[factor(samples$Treatment)], xlim=range(-2.5:2.5), ylim=range(-2.5:2.5), top=1000)

Group <- factor(paste(d$samples$group,samples$Age, sep="."))

cbind(d$samples, group=Group)

#GLM Test
design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)
my.contrasts <- makeContrasts(KOd4vsWTd4=KOd4-WTd4, KOd6vsWTd6=KOd6-WTd6, WTd4vsWTd6=WTd4-WTd6, KOd4vsKOd6=KOd4-KOd6, levels=design)
d <- estimateDisp(d, design)
fit <- glmQLFit(d, design)
qlf.KOd4.over.WTd4 <- glmQLFTest(fit, contrast=my.contrasts[,"KOd4vsWTd4"]);
qlf.KOd6.over.WTd6 <- glmQLFTest(fit, contrast=my.contrasts[,"KOd6vsWTd6"]);
qlf.WTd4.over.WTd6 <- glmQLFTest(fit, contrast=my.contrasts[,"WTd4vsWTd6"]);
qlf.KOd4.over.KOd6 <- glmQLFTest(fit, contrast=my.contrasts[,"KOd4vsKOd6"]);

qlf.KOd4.over.WTd4_tt <- topTags(qlf.KOd4.over.WTd4, Inf, sort.by="logFC");
qlf.KOd6.over.WTd6_tt <- topTags(qlf.KOd6.over.WTd6, Inf, sort.by="logFC");
qlf.WTd4.over.WTd6_tt <- topTags(qlf.WTd4.over.WTd6, Inf, sort.by="logFC");
qlf.KOd4.over.KOd6_tt <- topTags(qlf.KOd4.over.KOd6, Inf, sort.by="logFC");

write.table(qlf.KOd4.over.WTd4_tt, file = "qlf.KOd4.over.WTd4.txt", sep = "\t")
write.table(qlf.KOd6.over.WTd6_tt, file = "qlf.KOd6.over.WTd6.txt", sep = "\t")
write.table(qlf.KOd4.over.KOd6_tt, file = "qlf.KOd4.over.KOd6.txt", sep = "\t")
write.table(qlf.WTd4.over.WTd6_tt, file = "qlf.WTd4.over.WTd6.txt", sep = "\t")

qlf.KOd4.over.WTd4_tt_fdr <- qlf.KOd4.over.WTd4_tt[with(qlf.KOd4.over.WTd4_tt, order( qlf.KOd4.over.WTd4_tt$table$FDR)), ]
qlf.KOd6.over.WTd6_tt_fdr <- qlf.KOd6.over.WTd6_tt[with(qlf.KOd6.over.WTd6_tt, order( qlf.KOd6.over.WTd6_tt$table$FDR)), ]
qlf.KOd4.over.KOd6_tt_fdr <- qlf.KOd4.over.KOd6_tt[with(qlf.KOd4.over.KOd6_tt, order( qlf.KOd4.over.KOd6_tt$table$FDR)), ]
qlf.WTd4.over.WTd6_tt_fdr <- qlf.WTd4.over.WTd6_tt[with(qlf.WTd4.over.WTd6_tt, order( qlf.WTd4.over.WTd6_tt$table$FDR)), ]

write.table(qlf.KOd4.over.WTd4_tt_fdr, file = "qlf.KOd4.over.WTd4_GRCm39_tt_fdr.txt", sep = "\t")
write.table(qlf.KOd6.over.WTd6_tt_fdr, file = "qlf.KOd6.over.WTd6_GRCm39_tt_fdr.txt", sep = "\t")
write.table(qlf.KOd4.over.KOd6_tt_fdr, file = "qlf.KOd4.over.KOd6_GRCm39_tt_fdr.txt", sep = "\t")
write.table(qlf.WTd4.over.WTd6_tt_fdr, file = "qlf.WTd4.over.WTd6_GRCm39_tt_fdr.txt", sep = "\t")


rpkmTrue <- rpkm(d, normalized.lib.sizes= TRUE, log=TRUE);
write.table(rpkmTrue, file= "bea_2021_rpkm.txt", sep = "\t")

#write.table(rpmkTrue, file = "qlf.HL.over.HLC_rpkm.txt", sep = "\t")

cpmTRUE <- cpm(d, normalized.lib.sizes = TRUE, log=TRUE);
write.table(cpmTRUE, file= "bea_2021_cpm.txt", sep = "\t")



#Loads in trellis plot function
source("/Users/davehill/Desktop/RNASEQ/diffExpression_scripts/makeTrellisPlot.R")

#Loads in md plot function
source("/Users/davehill/Desktop/RNASEQ/diffExpression_scripts/makeMDPlot.R")

#Loads in volcano plot fuction
source("/Users/davehill/Desktop/RNASEQ/diffExpression_scripts/makeVolcanoPlot.R")

#Loads in make heatmap function
#takes a dge List object, a glm comparison object, "RPKM" or "CPM", a logical for Zscore setting, FDR level, clustering method, and linkage type
source("/Users/davehill/Desktop/RNASEQ/diffExpression_scripts/makeHeatMap.R")

#Loads in pca plot function
source("/Users/davehill/Desktop/RNASEQ/diffExpression_scripts/makePCAPlot.R")



#fuction call to make a trellis plot Inputs are d object and "CPM" or "RPKM"
#If type given is not CPM or RPKM then the function will automatically compute RPKM
makeTrellisPlot(d, "CPM")

#function call to make a MD plot
makeMDPlot(d, qlf.KOd4.over.WTd4, "KOd4", "WTd4")
makeMDPlot(d, qlf.KOd6.over.WTd6, "KOd6", "WTd6")
makeMDPlot(d, qlf.WTd4.over.WTd6, "WTd4", "WTd6")
makeMDPlot(d, qlf.KOd4.over.KOd6, "KOd4", "KOd6")

#function call to make a volcano plot
makeVolcanoPlot(qlf.KOd4.over.WTd4, "KOd4", "WTd4")
makeVolcanoPlot(qlf.KOd6.over.WTd6, "KOd6", "WTd6" )
makeVolcanoPlot(qlf.WTd4.over.WTd6, "WTd4", "WTd6")
makeVolcanoPlot(qlf.KOd4.over.KOd6, "KOd4", "KOd6")


#make enhanced volcano plot 
dataset = qlf.KOd4.over.WTd4_tt
dataset2 = qlf.KOd6.over.WTd6_tt
  
library(EnhancedVolcano)

EnhancedVolcano(dataset$table, lab = as.character(dataset$table$gene), x = 'logFC', y= 'PValue', title = dataset$comparison, selectLab = c("Il2ra", 'Cish','Ita','Il7r','Socs3','Tbx21'), boxedlabels = TRUE, drawConnectors =TRUE, widthConnectors = .75)


#function call to make a heatmap
makeHeatMap(counts, d, qlf.KOd4.over.WTd4, "RPKM", TRUE, .01)
makeHeatMap(counts, d, qlf.KOd6.over.WTd6, "CPM", TRUE, .35)
makeHeatMap(counts, d, qlf.WTd4.over.WTd6, "CPM", TRUE, .35)
makeHeatMap(counts, d, qlf.KOd4.over.KOd6, "CPM", TRUE, .35)



#function call to make a PCA plot
makePCAPlot(cpms)


#get intersection of up and down expressed genes
dataset = qlf.KOd4.over.KOd6_tt
dataset2 = qlf.WTd4.over.WTd6_tt
library(dplyr)

x = dataset$table[dataset$table$FDR < 0.05 & dataset$table$logFC > 1, ]
y = dataset2$table[dataset$table$FDR  < 0.05 & dataset$table$logFC <  -1, ]
both_genes_up <- intersect(x$gene, y$gene)
both_genes_down <- intersect(x$gene, y$gene)
KO_up_and_WT_down <-intersect(x$gene, y$gene)
KO_down_and_WT_up <-intersect(x$gene, y$gene)

