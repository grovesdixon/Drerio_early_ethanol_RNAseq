#wgcna4_module-correlations.R
#This script calculates and plots relationships between 
#module eigengenes and sample traits. Inputs come from:
#get_variance_stabilized_counts.R and wgcna3b_step-wise_network_construction.R
#code is adapted from examples given here: https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/


#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Load the expression and trait data saved in the first part
lnames = load(file = "wgcna/wgcna_input.Rdata");    
#The variable lnames contains the names of loaded variables.
lnames #(when you load a .Rdata object, if can have multiple variables assigned, the list of them is saved in this lnames variable, in this case the two objects you loaded are "datExpr" and "datTraits")


# Load network data saved in the second part.
lnames = load(file = "wgcna/wgcna3b_manual_sft15_minModSize10_cutDpth0.5_signed.Rdata") #This file using merged dynamic modules. The ethanol modules are largely retained, but the two that are downregulated are merged into one.
lnames    #the outputs from network construction: "MEs"          "moduleLabels" "moduleColors" "geneTree"

#match up the rownames in the expression 
sampleNames = datTraits$sample.names
rownames(datTraits) = sampleNames  #make the sample names the rownames for dataframe datTraits
datTraits$sample.names<-NULL #get rid of the sample.names column
rownames(datExpr)     #this tells us the rownames as given to WGCNA
rownames(datTraits)   #this is for the entire set of datTraits and includes samples that were not put into WGCNA
rownames(MEs) = rownames(datExpr)
datTraits = datTraits[rownames(datTraits) %in% rownames(datExpr),]  #reduce datTraits to those we have WGCNA results for
#if everything is matched up this should say TRUE
sum(rownames(datTraits) == rownames(datExpr)) == length(rownames(datExpr))
colnames(datTraits) = c('ethanol', 'time', 'seqjob', 'rand.eth')
head(datTraits)
dim(datTraits)

datTraits$ethanol = as.numeric(datTraits$ethanol)
datTraits$rand.eth = as.numeric(datTraits$rand.eth)



#add in binaries for timepoints
datTraits = datTraits %>% 
  mutate(time=as.numeric(as.character(time)),
         seqjob=as.numeric(seqjob),
         `8hpf` = if_else(time==8,
                         1,
                         0),
         `10hpf` = if_else(time==10,
                       1,
                       0),
         `14hpf` = if_else(time==14,
                            1,
                            0)) %>% 
  dplyr::select(ethanol, time, `8hpf`, `10hpf`, `14hpf`, seqjob, rand.eth)
rownames(datTraits)=sampleNames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate module eigengenes with color labels
#This will provide the first principal component for expression
#behavior for the genes in each module. On that principal component,
#each sample will have a loading value. These values can be corelated
#with our known sample traits to get an idea what biological mechanism
#the co-reguated genes in a given module might be responding to
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#use the cor() function to get the correlations between the module eigengenes and the trait data
moduleTraitCor = cor(MEs, datTraits, use = "p");
#get p values as well
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#write out module loadings
mEigs=MEs
rownames(mEigs) = rownames(datTraits)
# save(mEigs, file='wgcna/moduleEigengenes.Rdata')

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)

#now replot the heatmap
sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(subCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)

#set up additional formatting variables
rows = rownames(moduleTraitCor)
sub.colors = substr(rownames(moduleTraitCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = rownames(moduleTraitCor),
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#compare p-values for read and randomized ethanol
plot(density(moduleTraitPvalue[,'ethanol']), lwd=2)
lines(density(moduleTraitPvalue[,'rand.eth']), lwd=2, col='red')




#estimate false discovery rate
BREAKS=40
hist(moduleTraitPvalue[,'ethanol'], breaks= BREAKS, col='blue', main='', xlab='P-value')
hist(moduleTraitPvalue[,'rand.eth'], breaks= BREAKS, add=T, col='red')
legend('topright', c('real', 'false'), pt.bg=c('blue', 'red'), pch=22)
abline(v=0.02, lwd=2,lty=2)


# get time hub genes ------------------------------------------------------
#get all modules significantly associated with time, or any particular time point
library(tidyverse)
CUT=0.05
timeMods = moduleTraitPvalue %>% 
  as_tibble() %>% 
  mutate(module=rownames(moduleTraitPvalue)) %>% 
  filter(time<CUT)
dim(timeMods)

#get module membership
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")) %>% 
  select(rownames(moduleTraitCor))
head(geneModuleMembership)

#save both
save(geneModuleMembership, timeMods, file='results/timeMods.Rdata')


#-------- replot the heatmap after subsetting for modules that are significant for your trait of interest--------
#select a significance cutoff (CUT) and the trait you are interested in (TRAIT)
CUT = 0.05
TRAIT = 'ethanol'
# TRAIT = 'time'
# TRAIT = 'seqjob'


#subset the correlations
subCor = moduleTraitCor[moduleTraitPvalue[, TRAIT] < CUT,]
subP = moduleTraitPvalue[moduleTraitPvalue[, TRAIT] < CUT,]
rows = rownames(moduleTraitCor)[moduleTraitPvalue[, TRAIT] < CUT]
sub.colors = substr(rownames(subCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")

#now replot the heatmap
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(subCor, 2), "\n(",
                           signif(subP, 1), ")", sep = "");
dim(textMatrix) = dim(subCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
COLORS = greenWhiteRed(50)
blueWhiteRed(50)


labeledHeatmap(Matrix = subCor,
               xLabels = names(datTraits),
               yLabels = rownames(subCor),
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
               
               
#=====================================================================================
#REPLOT THE MODULE EIGENGENE CLUSTERING TREE
#This shows you how the modules are similar to one another.
#Pushing the merging theshold (argument 3 for wgcna3b_step-wise_network_construction.R) will join modules together
#This is like moving a horizontal line up this tree figure, if the line is above a node, modules below that node will be joined into one
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

#plot them
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================

#=====================================================================================
#
#  Code chunk 4 Gather module membership data for the genes in each module
#
#=====================================================================================
#"module membership" is a measure of how strongly individual genes correlate with the 
#module eigengene. The gene that matches best with the module eigengene can be thought
#of as a hub gene, (ie it's variation in expression across the samples is most exemplary 
#of the module)

# Define dataframe trait.df containing the a trait of interest from datTraits
trait.df = as.data.frame(datTraits[,TRAIT]);
names(trait.df) = TRAIT
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
head(geneModuleMembership)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, trait.df, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(trait.df), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.df), sep="");
modules = sub.colors
plot.cols = modules 

#=====================================================================================
#
#  Code chunk 5 plot scatterplots of module membership and trait correlations
#
#=====================================================================================   
#The code below loops through each of the modules that were significant for the TRAIT of 
#interest that you chose above. Each point in the scatterpot is a gene in the module. 
#For each one of the it plots the the genes' correlation
#with the trait of interest against the genes' module memberships,
#(the correlation between the gene's variation accross samples and the module eigengene)
#A tight correlation here is suggestive that the correlated variation in gene expression 
#captured by the module is truely associated with the trait of interest.
#A tight correlation basically means, the better a gene fits into this module, the more strongly
#it correlates with the trait.



quartz()
par(mar=c(5,5))
m='darkolivegreen4'
for (m in modules){
	column = match(m, modNames);
	moduleGenes = moduleColors==m;
	
	# sizeGrWindow(7, 7);
	# par(mfrow = c(1,1));
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                   abs(geneTraitSignificance[moduleGenes, 1]),
	                   xlab = paste("Module Membership in", m, "module"),
	                   ylab = paste("Correlation with", TRAIT),
	                   main = paste("Module membership vs. gene significance\n"),
	                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black', bg = m, pch = 21, cex = 1.5)
}


###################################################
#### LOOK AT RESULTS FOR A A PARTICULAR MODULE ####
###################################################

source("deseq/zebrafish_RNAseq_functions.R")
library("biomaRt")
embl = useMart("ensembl", dataset="drerio_gene_ensembl")
# listAttributes(embl)


##### BUILD BOXPLOTS #####
#gather genes for modules of interest
genes = rownames(geneModuleMembership)
modules = c('mediumpurple4', 'darkolivegreen4')
moduleGenes = c()
mod = c()
for (m in modules){
	print("---------------------")
	print(m)
	toAdd=genes[moduleColors==m]
	mlabs=rep(m, length(toAdd))
	print(paste(length(toAdd), "genes added"))
	print(toAdd)
	print(mlabs)
	moduleGenes = append(moduleGenes, toAdd)
	mod=append(mod, mlabs)
}
mres=data.frame(moduleGenes, mod)
rownames(mres) = mres[,1]
head(mres)

#merge with DESeq results
lnames=load("deseq/ethanol_full_LRT_results.Rdata")
head(res.eth)
r=data.frame(res.eth[rownames(res.eth) %in% moduleGenes,])
mdat=merge(mres,r,by=0)
mdat
nrow(mdat) == nrow(mres)

# b=boxplot(mdat$log2FoldChange~mdat$mod, axes=F, notch=T, outline=F, ylab="Log2 Fold Difference", ylim=c(-2,2.1))
mdat$xs=match(mdat$mod, rev(modules))
plot(mdat$log2FoldChange~jitter(mdat$xs, factor=1.2), bg=mdat$mod, pch=21, cex=1.2, ylim=c(-2,2.1), xlim=c(.5,2.5), axes=F, ylab="Log2 Fold Difference", xlab='')
axis(2)
axis(1, at=1:length(modules), label=rev(modules))
b=boxplot(mdat$log2FoldChange~mdat$mod, axes=F, outline=F, add=T, lwd=1.5, xlab="")
abline(h=0, lty=2, col='black', lwd=1.5)


#do ggplot boxplots
library(ggplot2)
p <- ggplot(mdat, aes(x=mod, y= log2FoldChange)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) + theme_bw() + xlab("Module") + ylab("Log2 Fold Difference")
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, fill=rev(mod)) + geom_hline(yintercept=0, linetype='dashed', color='black', size=1)




#pick the module to look at
m='mediumpurple4'
m='darkolivegreen4'



#gather gene names for this module
genes = rownames(geneModuleMembership)
moduleGenes = genes[moduleColors==m]
moddat = data.frame(geneModuleMembership[moduleGenes, paste("MM",m,sep='')]);rownames(moddat) = moduleGenes;colnames(moddat)=c(m)
head(moddat)
dat=merge_gene_names(moddat, sort.column=m)
head(dat)
length(moduleGenes)
dim(dat)
dat


#write out the gene names
fileName=paste(m, 'gene_data.tsv',sep="_")
write.table(dat, file=paste('results', fileName, sep="/"), sep="\t", quote=F, row.names=F)


#LOOK AT DISTRIBUTION ACCROSS CHROMOSOME 4
datc4=na.omit(dat[dat$chromosome_name==4,])
print("Percentage of module on Chromosome 4:")
nrow(datc4)/nrow(dat)*100
chrom4_length = 76625712 #taken from here http://www.ensembl.org/Danio_rerio/Location/Chromosome?r=4:59726329-59826329
plot(density(datc4$start_position/1e6, bw=4), xlim=c(0,chrom4_length/1e6+1), main="Module Gene Locations", xlab='Chromosome 4 position (Mb)', axes=F)
axis(1, at = seq(10, 70, by = 10))
axis(2)


#BUILD HEATMAPS TO SHOW EXPRESSION DIFFERENCES
lnames=load("./datasets/large_ignored/raw_rld.Rdata")
head(rld.df)

#select for this module
x=rld.df[rownames(rld.df) %in% moduleGenes, colnames(rld.df) %in% sampleNames]
dim(x)
length(moduleGenes)
head(x)



#sort the dataframe by ethanol treatment, and time
datTraits$ordTime = as.numeric(as.character(datTraits$time))
y=datTraits[with(datTraits, order(ordTime, ethanol)), ]
y=datTraits[with(datTraits, order(ethanol, ordTime)), ]
# y=y[y$time==2,]
x=x[,rownames(y)]
colnames(x)


#get z-scores for genes
head(x)
gmeans = apply(x, 1, mean)
gsds = apply(x, 1, sd)

z=(x-gmeans)/gsds
head(z)

library(pheatmap)
library(plotrix)
labs = get_gene_names(rownames(z))
labels = labs$description
for (i in 1:length(labels)){
	if(labels[i] == ""){
		labels[i]<-labs$external_gene_name[i]
	}
}

pheatmap(z,cluster_cols=F, treeheight_row=0, border_color=NA,clustering_distance_rows="correlation", labels_row=labels)





