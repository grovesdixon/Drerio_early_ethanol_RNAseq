#initialize_counts_step2.R

library(sva)
library(ape)
library(vegan)

ll=load("datasets/large_ignored/initialize_countsImage1.Rdata")

# IDENTIFY POORLY SEQUENCED GENES AND OUTLIERS AS IN WGCNA TUTORIAL -------


#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================

datExpr0 = as.data.frame(t(rld.df));

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have genes or samples with too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
head(gsg)


#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
rld.df=t(datExpr0)
rld=rld[rownames(rld.df),]
dim(rld.df)
dim(rld)
nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#build sample heatmaps to confirm the outlier results
library(pheatmap)
quartz()
pheatmap(cor(rld.df))
x=colnames(rld.df)
y=sub("NoE", "C", x)
z=sub("h.", "", y)
l = substr(z, 1,nchar(z)-1)
l2=sub(".", "", l, fixed=T)
# labs=paste(l2, seqjob, sep="_") #to include batch
labs=(l2)
pheatmap(cor(rld.df, method = 'spearman'), labels_row=labs, labels_col=labs, treeheight_row=0, treeheight_col=0, number_color='blue')


#=====================================================================================
#
#  Code chunk 6
# 
#=====================================================================================

#Remove outliers by setting a branch cut threshold
# Plot a line to show the cut
cut.height = 150
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cut.height, col = "red", lty = 2);
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
keepSampleNames = rownames(datExpr0)[keepSamples]
outlierNames = rownames(datExpr0)[clust==0]
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #number of samples left after outlier removal
print(paste(length(outlierNames), "samples were flagged as outliers and removed:"))
outlierNames
print(paste(nSamples, "samples were kept"))

#save the outlier names so you can optionally remove them in other analyses
# save(outlierNames, file = 'datasets/outliers.Rdata')

#=====================================================================================



# MULTIVARIATE ANALYSES  --------------------------------------------------


#here we want to explore the dataset to decide what to do about batch effects and different treatments
source("./initialize_counts/multivariate_functions.R")

#PCA WITH DESEQ FUNCTION
library(ggplot2)
library(DESeq2)
NTOP = 25000
point.size = 5
plotPCA(rld,intgroup=c("treatment"),ntop= NTOP) #run with original DESeq function (just to prove they work the same)
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment") #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'time', main = 'Time (hr)')
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'seqjob', main="Batch")


#look at PCs 3 -10:
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=3, pc2=4) #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=5, pc2=6) #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=7, pc2=8) #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=9, pc2=10) #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=1, pc2=8) #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=8, pc2=9) #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment", pc1=1, pc2=9) #with modified one

#BUILD THE SAME PLOTS AFTER REMOVING OUTLIERS FLAGGED ABOVE
#(this is a remnant from when we had outliers, don't anymore so it changes nothing)
#remove outliers from rld.df (not that rld.df and datExpr are just transpositions of one another)
dim(rld.df)
rld.df = rld.df[, !colnames(rld.df) %in% outlierNames]
dim(rld.df) #(removed 6 outlier samples, from 55 down to 49)
coldata = coldata[!coldata$sample.names %in% outlierNames,]

#replot
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment") #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'time', main = 'Time (hr)')
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'seqjob', main="Batch")
#note that we get much better clustering by time once outliers were removed
#supports that these really were just badly sequenced

#---------------- REMOVE 6 HR SAMPLES ----------------#
#these are isolated within seqjob1 and do not have ethanol treated replicates
#hence they likely add a lot of variation to dataset we are not interested in
#and will confuse the batch adjustment
#remove them before filtering for variation and for batch effects
head(rld.df)
dim(rld.df)
sixHourSamples=colnames(counts)[grep('6h', colnames(counts))]
length(sixHourSamples)
rld.df=rld.df[, !colnames(rld.df) %in% sixHourSamples]     # (remove 4 six hour samples (one was an outlier))
coldata=coldata[!coldata$sample.names %in% sixHourSamples,]# 45 out of original of 55 samples remaining
dim(rld.df)
dim(coldata)

#---------------- BASE MEAN FILTERING ----------------#
counts.no.out = counts[,!colnames(counts) %in% outlierNames]
counts.no.out = counts.no.out[,!colnames(counts.no.out) %in% sixHourSamples]
dim(counts.no.out)
dim(rld.df)


CUT=5
table(res$baseMean>CUT)
keep=rownames(res)[res$baseMean>CUT]
length(keep)
rld.df=rld.df[keep,]
dim(rld.df)

#replot
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'treatment', main = "Ethanol Treatment") #with modified one
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'time', main = 'Time (hr)')
mod.plotPCA.df(df = rld.df, coldat = coldata, intgroup = 'seqjob', main="Batch")


#USE VEGAN TO PARTITION VARIANCE
library(ape)
library(vegan)
datExpr=t(rld.df)
coldata$randomized.ethanol=sample(coldata$treatment, replace=F)
sum(coldata$sample.names==rownames(datExpr)) == nrow(datExpr) #double-check coldata and datExpr match
ad=adonis(datExpr~seqjob+time+treatment+ randomized.ethanol,data=coldata,method="manhattan")
labs=c("batch","time", "ethanol", "random ethanol", "residuals")
cols=c("skyblue","green2","coral", "black", "grey80")
pie(ad$aov.tab$R2[1:5],labels=labs,col=cols,main="Variance Partitioning")
print(ad)


#SUBSET DATA TO LOOK FOR ETHANOL EFFECT
#look at controls alone
samples = colnames(rld.df)
controls = append(samples[grep('C', samples)], samples[grep('NoE', samples)])
control.df = rld.df[, colnames(rld.df) %in% controls]
control.coldat = coldata[coldata$sample.names %in% controls,]
eth.df = rld.df[, !colnames(rld.df) %in% controls]
eth.coldat = coldata[!coldata$sample.names %in% controls,]
mod.plotPCA.df(df = control.df, coldat = control.coldat, intgroup = 'time', main="Time (Control samples only)")
mod.plotPCA.df(df = eth.df, coldat = eth.coldat, intgroup = 'time', main="Time (Ethanol samples only)")

#look within sequencing jobs (batches)
sum(colnames(rld.df) == coldata$sample.names)==ncol(rld.df)
j1=rld.df[, colnames(rld.df) %in% coldata$sample.names[coldata$seqjob==1]]
j2=rld.df[, colnames(rld.df) %in% coldata$sample.names[coldata$seqjob==2]]
j1c=coldata[coldata$seqjob==1,]
j2c=coldata[coldata$seqjob==2,]
mod.plotPCA.df(df = j1, coldat = j1c, intgroup = 'treatment', main="Treatment (Batch1 only)", pc1=1, pc2=2)
mod.plotPCA.df(df = j2, coldat = j2c, intgroup = 'treatment', main="Treatment (Batch2 only)", pc1=1, pc2=2)
mod.plotPCA.df(df = j1, coldat = j1c, intgroup = 'time', main="Time (Batch1 only)")
mod.plotPCA.df(df = j2, coldat = j2c, intgroup = 'time', main="Time (Batch2 only)")

#----------- LOOK WITHIN TIMEPOINTS -----------#
#8 HOUR SAMPLES
timepoint=8
eight=rld.df[,colnames(rld.df) %in% coldata$sample.names[coldata$time==timepoint]]
ec=coldata[coldata$time==timepoint,]
#overall treatment and batch effects within 8-hour timepoint accross both experiments
mod.plotPCA.df(df = eight, coldat = ec, intgroup = 'treatment', main=bquote("Treatment ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)
mod.plotPCA.df(df = eight, coldat = ec, intgroup = 'seqjob', main=bquote("Batch ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)

#adonis
ec$randomized.ethanol=sample(ec$treatment, replace=F)
ad=adonis(t(eight)~seqjob+treatment+randomized.ethanol,data=ec,method="manhattan")
labs=c("batch", "ethanol", "random ethanol", "residuals")
cols=c("skyblue","green2", "black", "grey80")
pie(ad$aov.tab$R2[1:4],labels=labs,col=cols,main="Variance Partitioning")
print(ad)

#subset for batches
e1=eight[, colnames(eight) %in% ec$sample.names[ec$seqjob==1]]
e2= eight[, colnames(eight) %in% ec$sample.names[ec$seqjob==2]]
e1c=ec[ec$seqjob==1,]
e2c=ec[ec$seqjob==2,]
#look at ethanol effect within batches within timepoint 8
mod.plotPCA.df(df = e1, coldat = e1c, intgroup = 'treatment', main=bquote("Treatment (job1 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)
mod.plotPCA.df(df = e2, coldat = e2c, intgroup = 'treatment', main=bquote("Treatment (job2 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)
#adonis
e2c$randomized.ethanol=sample(e2c$treatment)
ad=adonis(t(e2)~treatment+randomized.ethanol,data=e2c,method="manhattan")
labs=c("ethanol", "random ethanol", "residuals")
cols=c("green2", "black", "grey80")
pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="Variance Partitioning")
print(ad)
#only marginal effect apparent


#10 HOUR SAMPLES
timepoint=10
ten=rld.df[,colnames(rld.df) %in% coldata$sample.names[coldata$time==timepoint]]
ec=coldata[coldata$time==timepoint,]
#overall batch and ethanol in timepoint 10
mod.plotPCA.df(df = ten, coldat = ec, intgroup = 'treatment', main=bquote("Treatment ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)
mod.plotPCA.df(df = ten, coldat = ec, intgroup = 'seqjob', main=bquote("Batch ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)

#adonis
ec$randomized.ethanol=sample(ec$treatment, replace=F)
ad=adonis(t(ten)~seqjob+treatment+randomized.ethanol,data=ec,method="manhattan")
labs=c("batch", "ethanol", "random ethanol", "residuals")
cols=c("skyblue","green2", "black", "grey80")
pie(ad$aov.tab$R2[1:4],labels=labs,col=cols,main="Variance Partitioning")
print(ad)

#subset for batches
e1=ten[, colnames(ten) %in% ec$sample.names[ec$seqjob==1]]
e2= ten[, colnames(ten) %in% ec$sample.names[ec$seqjob==2]]
e1c=ec[ec$seqjob==1,]
e2c=ec[ec$seqjob==2,]
e2c$randomized.ethanol=sample(e2c$treatment)
# mod.plotPCA.df(df = e1, coldat = e1c, intgroup = 'treatment', main=bquote("Treatment (job1 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)#not enough samples because of outlier removal
mod.plotPCA.df(df = e1, coldat = e1c, intgroup = 'treatment', main=bquote("Treatment (job1 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)
mod.plotPCA.df(df = e2, coldat = e2c, intgroup = 'treatment', main=bquote("Treatment (job2 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)

#adonis
ad=adonis(t(e2)~treatment+randomized.ethanol,data=e2c,method="manhattan")
labs=c("ethanol", "random ethanol", "residuals")
cols=c("green2", "black", "grey80")
pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="Variance Partitioning")
print(ad)
#significant effect of treatment within timepoint 10 batch 2, but not batch1


#14 hour samples
timepoint=14
forteen=rld.df[,colnames(rld.df) %in% coldata$sample.names[coldata$time==timepoint]]
ec=coldata[coldata$time==timepoint,]
mod.plotPCA.df(df = forteen, coldat = ec, intgroup = 'treatment', main=bquote("Treatment ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)
mod.plotPCA.df(df = forteen, coldat = ec, intgroup = 'seqjob', main=bquote("Batch ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)
e1=forteen[, colnames(forteen) %in% ec$sample.names[ec$seqjob==1]]
e2= forteen[, colnames(forteen) %in% ec$sample.names[ec$seqjob==2]]
e1c=ec[ec$seqjob==1,]
e2c=ec[ec$seqjob==2,]
e2c$randomized.ethanol=sample(e2c$treatment)
mod.plotPCA.df(df = e2, coldat = e2c, intgroup = 'treatment', main=bquote("Treatment (job2 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)

#adonis
ec$randomized.ethanol=sample(ec$treatment, replace=F)
ad=adonis(t(forteen)~treatment+randomized.ethanol,data=ec,method="manhattan")
labs=c("ethanol", "random ethanol", "residuals")
cols=c("green2", "black", "grey80")
pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="Variance Partitioning")
print(ad)


#end of data exploration.
#conclusions:
# 1. lots of unexplained variation
# 2. time is significant and strongest effect
# 3. batch is significant and second strongest
# 4. ethanol is at best marginally significant when time and batch have been accounted for

#=========================================================================


# CONTROLLING FOR BATCH EFFECTS -------------------------------------------


#adjust for batch effects using combat in sva package
sum(coldata$sample.names == colnames(rld.df)) == ncol(rld.df)

#adonis to see partitioning before adjustment
coldata$randomized.ethanol=sample(coldata$treatment, replace=F)
ad=adonis(t(rld.df)~seqjob+time+treatment+randomized.ethanol,data=coldata,method="manhattan")
labs=c("batch","time", "ethanol", "random ethanol", "residuals")
cols=c("skyblue","green2","coral", "black", "grey80")
pie(ad$aov.tab$R2[1:5],labels=labs,col=cols,main="Variance Partitioning")
print(ad)

#batch is marginally significant, and we saw in PCAs above that
#we cannot see separation of ethanol samples until we look inside of batches
#based on this will adjust for batch effects before outputting for WGCNA


#assign the batch variable
batch=coldata$seqjob

#create a model matrix for the adjustment variables, including the variable of interest. Note that you do not include batch in creating this model matrix - it will be included later in the ComBat function. In this case there are no other adjustment variables so we simply fit an intercept term.
modcombat = model.matrix(~1, data=coldata)

# Note that adjustment variables will be treated as given to the ComBat function. This means if you are trying to adjust for a categorical variable with p different levels, you will need to give ComBat p-1 indicator variables for this covariate. We recommend using the model.matrix function to set these up. For continuous adjustment variables, just give a vector in the containing the covariate values in a single column of the model matrix.
# We now apply the ComBat function to the data, using parametric empirical Bayesian adjustments.
cbat.rld = ComBat(dat=rld.df, batch=batch, mod=modcombat, par.prior=TRUE)
head(cbat.rld)
dim(rld.df)
dim(cbat.rld)

#adonis again
coldata$randomized.ethanol=sample(coldata$treatment, replace=F)
ad=adonis(t(cbat.rld)~seqjob+time+treatment+randomized.ethanol,data=coldata,method="manhattan")
labs=c("batch","time", "ethanol", "random ethanol", "residuals")
cols=c("skyblue","green2","coral", "black", "grey80")
pie(ad$aov.tab$R2[1:5],labels=labs,col=cols,main="Variance Partitioning")
print(ad)
#time significant, batch effect removed, ethanol marginally significant, ranodom ethanol usually ns

#------------ SAVE RESULTS FOR WGCNA ------------
#change name to fit with WGCNA tutorial
datTraits = coldata
datExpr = t(cbat.rld)
dim(datExpr)
dim(datTraits)

# datExpr=t(rld.df)

#save data for full dataset
save(datExpr, datTraits, file = "wgcna/wgcna_input.RData") #use the Rdata file for input to run WGCNA on full dataset
