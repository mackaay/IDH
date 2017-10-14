library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# Read the data into R
seqdata <- read.delim("data.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("sampleinfo.txt")
head(seqdata)
dim(seqdata)


countdata <- seqdata
# Look at the output
head(countdata)
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
countdata <- as.matrix(countdata)

barplot(countdata,names=colnames(countdata))
# Add a title to the plot
title("Barplot of library sizes")
logcounts <- countdata
logcounts <- as.matrix(logcounts)
class(logcounts) <- "numeric"
logcounts <- logcounts[,-1]
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts")
abline(h=median(logcounts),col="red")
title("Boxplots of miRNAs across IDH mutant and wild type GBM")

#ggfortify#
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
install_github('sinhrks/ggfortify')
library(ggfortify)


pca1<- prcomp(clinical_data_to_TCGA_miRNA_IDH1_mutation_master[,-c(1,2)])
as.factor(clinical_data_to_TCGA_miRNA_IDH1_mutation_master$IDH1status)
autoplot(pca1, data = clinical_data_to_TCGA_miRNA_IDH1_mutation_master, 
         colour = 'IDH1status', size = 3) + scale_color_brewer(palette='Set1')+ theme_bw()

# How many cell types and in what order are they stored?
levels(sampleinfo$IDH1mutation)
# Let's choose purple for basal and orange for luminal
col.idh <- c("red","blue", "green", "purple", "orange")[sampleinfo$IDH1mutation]
data.frame(sampleinfo$IDH1mutation,col.idh)


##Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:50]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable genes across samples",ColSideColors=col.idh,scale="row")
# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 50 most variable genes across samples",ColSideColors=col.idh,scale="row")
dev.off()



y <- DGEList(logcounts)
plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$IDH1mutation)
plotMDS(y,col=col.idh)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group <- paste(sampleinfo2$IDH1mutation)
group <- factor(group)
##group matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
rownames(design) <- colnames(logcounts)
design

##contrast matrix
contrast.matrix <- makeContrasts(paste0(unique(group), collapse = "-"), levels = design)
contrast.matrix

##step1
fit <- lmFit(logcounts, design)
##step2
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)
dim(fit2)
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
idhDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(idhDEG)
write.csv(idhDEG, file = "idhDEmiR.csv")




# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit2,coef=1,highlight=100)


par(mfrow=c(1,1))
with(idhDEG, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", xlim=c(-3,3)), cex = .9 )
with(subset(idhDEG, abs(logFC)>1), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
library(calibrate)
idhDEG$name <- rownames(idhDEG)
with(subset(idhDEG, abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs=name, cex=.9))


idhDEG_select <- c("hsa-miR-204", "hsa-miR-34b", "hsa-miR-221","hsa-miR-148a", 
                   "hsa-miR-155", "hsa-miR-34a", "hsa-miR-222")
idhDEG_select
idhDEG_8miR <- logcounts[idhDEG_select,]
dim(idhDEG_8miR)
head(idhDEG_8miR)

heatmap.2(idhDEG_8miR,col=rev(morecols(50)),trace="none", 
          main="IDH DE miRs",ColSideColors=col.idh,scale="row", margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$IDH1mutation), col = unique(col.idh), lty = 1, lwd= 5, cex=.7)



