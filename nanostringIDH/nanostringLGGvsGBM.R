library(edgeR)
library(limma)
library(gplots)
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
# Store gene id as rownames
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
title("Boxplots of miRNAs over all serum cohort samples")

#ggfortify#
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
install_github('sinhrks/ggfortify')
library(ggfortify)

###loading data#####
pca1<- prcomp(RUVIII_normalized_data_transformed[,-c(1,2,3,4,5,6)])
as.factor(RUVIII_normalized_data_transformed$Grade)
autoplot(pca1, data = RUVIII_normalized_data_transformed, 
         colour = 'Grade', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()


levels(sampleinfo$Grade)
# Let's choose purple for basal and orange for luminal
col.grade <- c("blue", "green")[sampleinfo$Grade]
data.frame(sampleinfo$Grade,col.grade)


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
          main="Top 50 most variable genes across samples",ColSideColors=col.grade,scale="row",
          margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$Grade), col = unique(col.grade), lty = 1, lwd= 5, cex=.7)



# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable genes across samples",ColSideColors=col.idh,scale="row",
          margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$IDHstatus), col = unique(col.idh), lty = 1, lwd= 5, cex=.7)
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
group <- paste(sampleinfo$Grade)
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
gradeDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(gradeDEG)
write.csv(gradeDEG, file = "gradeDEmiR.csv")

##################


# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit2,coef=1,highlight=100)


par(mfrow=c(1,1))
with(gradeDEG, plot(logFC, -log(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2.5)), cex = .9 )
with(subset(gradeDEG, abs(logFC)>1), points(logFC,-log(P.Value), pch=20, col="red"))
library(calibrate)
gradeDEG$name <- rownames(gradeDEG)
with(subset(gradeDEG, abs(logFC)>1),  textxy(logFC, -log(P.Value), labs=name, cex=.9))


gradeDEG_select <- c("hsa-miR-451a", "hsa-miR-223-3p", "hsa-miR-548aa+hsa-miR-548t-3p")
gradeDEG_select
gradeDEG_3miR <- logcounts[gradeDEG_select,]
dim(gradeDEG_3miR)
head(gradeDEG_3miR)

heatmap.2(gradeDEG_3miR,col=rev(morecols(50)),trace="none", 
          main="Grade DE miRs",ColSideColors=col.grade,scale="row", 
          margins = c(5,15), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$Grade), col = unique(col.grade), lty = 1, lwd= 5, cex=.7)


######batch effect??#####
levels(sampleinfo$batch)
col.batch <- c(brewer.pal(8,"Accent"))[sampleinfo$batch]
col.batch

heatmap.2(gradeDEG_3miR,col=rev(morecols(50)),trace="none", 
          main="Grade DE miRs",ColSideColors=col.batch,scale="row", 
          margins = c(5,15), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$batch), col = unique(col.batch), lty = 1, lwd= 5, cex=.7)
####################






###################################Remove batch 3 and 4#############################
seqdata2 <- read.delim("data2.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo2 <- read.delim("sampleinfo2.txt")
head(seqdata2)
dim(seqdata2)


countdata2 <- seqdata2
# Look at the output
head(countdata2)
# Store gene id as rownames
rownames(countdata2) <- seqdata2[,1]
head(countdata2)
colnames(countdata2)
countdata2 <- as.matrix(countdata2)

barplot(countdata2,names=colnames(countdata2))
# Add a title to the plot
title("Barplot of library sizes")
logcounts2 <- countdata2
logcounts2 <- as.matrix(logcounts2)
class(logcounts2) <- "numeric"
logcounts2 <- logcounts2[,-1]
# Check distributions of samples using boxplots
boxplot(logcounts2, xlab="", ylab="Log2 counts")
abline(h=median(logcounts),col="red")
title("Boxplots of miRNAs over 74 serum cohort samples")


#ggfortify#
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
install_github('sinhrks/ggfortify')
library(ggfortify)

###loading data#####
pca1<- prcomp(RUVIII_normalized_data_transformed2[,-c(1,2,3,4,5,6,7)])
as.factor(RUVIII_normalized_data_transformed2$Grade)
autoplot(pca1, data = RUVIII_normalized_data_transformed2, 
         colour = 'Grade', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()
as.factor(RUVIII_normalized_data_transformed2$IDH1)
autoplot(pca1, data = RUVIII_normalized_data_transformed2, 
         colour = 'IDH1', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()
as.factor(RUVIII_normalized_data_transformed2$batch)
autoplot(pca1, data = RUVIII_normalized_data_transformed2, 
         colour = 'batch', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()



levels(sampleinfo2$Grade)
# Let's choose purple for basal and orange for luminal
col.grade2 <- c("blue", "green")[sampleinfo2$Grade]
data.frame(sampleinfo2$Grade,col.grade2)


##Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes2 <- apply(logcounts2, 1, var)
head(var_genes2)
dim(var_genes2)
# Get the gene names for the top 500 most variable genes
select_var2 <- names(sort(var_genes2, decreasing=TRUE))[1:50]
head(select_var2)
# Subset logcounts matrix
highly_variable_lcpm2 <- logcounts2[select_var2,]
dim(highly_variable_lcpm2)
head(highly_variable_lcpm2)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
heatmap.2(highly_variable_lcpm2,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable genes across samples",ColSideColors=col.grade2,scale="row",
          margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo2$Grade), col = unique(col.grade2), lty = 1, lwd= 5, cex=.7)



# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable genes across samples",ColSideColors=col.idh,scale="row",
          margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$IDHstatus), col = unique(col.idh), lty = 1, lwd= 5, cex=.7)
dev.off()



y <- DGEList(logcounts2)
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
group2 <- paste(sampleinfo2$Grade)
group2 <- factor(group2)
##group matrix
design2 <- model.matrix(~0+group2)
colnames(design2) <- levels(group2)
rownames(design2) <- colnames(logcounts2)
design2

##contrast matrix
contrast.matrix2 <- makeContrasts(paste0(unique(group2), collapse = "-"), levels = design2)
contrast.matrix2

##step1
fit_2 <- lmFit(logcounts2, design2)
##step2
fit2_2<- contrasts.fit(fit_2, contrast.matrix2)
fit2_2<- eBayes(fit2_2)
dim(fit2_2)
##step3
tempOutput2 = topTable(fit2_2, coef=1, n=Inf)
gradeDEG2 = na.omit(tempOutput2) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(gradeDEG2)
write.csv(gradeDEG2, file = "gradeDEmiR2.csv")

##################


# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit2,coef=1,highlight=100)


par(mfrow=c(1,1))
with(gradeDEG2, plot(logFC, -log(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2.5)), cex = .9 )
with(subset(gradeDEG2, abs(logFC)>1), points(logFC,-log(P.Value), pch=20, col="red"))
library(calibrate)
gradeDEG2$name <- rownames(gradeDEG2)
with(subset(gradeDEG2, abs(logFC)>1),  textxy(logFC, -log(P.Value), labs=name, cex=.9))


gradeDEG_select2 <- c("hsa-miR-451a", "hsa-miR-223-3p")
gradeDEG_select2
gradeDEG_2miR <- logcounts2[gradeDEG_select2,]
dim(gradeDEG_2miR)
head(gradeDEG_2miR)

heatmap.2(gradeDEG_2miR,col=rev(morecols(50)),trace="none", 
          main="Grade DE miRs",ColSideColors=col.grade2,scale="row", 
          margins = c(5,15), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo2$Grade), col = unique(col.grade), lty = 1, lwd= 5, cex=.7)



######batch effect??#####
levels(sampleinfo2$batch)
col.batch2 <- c(brewer.pal(7,"Accent"))[sampleinfo2$batch]
col.batch2

heatmap.2(gradeDEG_2miR,col=rev(morecols(50)),trace="none", 
          main="Grade DE miRs",ColSideColors=col.batch2,scale="row", 
          margins = c(5,15), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$batch), col = unique(col.batch), lty = 1, lwd= 5, cex=.7)
####################





