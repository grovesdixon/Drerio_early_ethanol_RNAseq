#DESeq_early_ethanol.R
#Groves Dixon
#written 8-18-16
#updated 4-18-19
#comments in quotes taken from here: http://www.bioconductor.org/help/workflows/rnaseqGene/#time-course-experiments
#use this script to test for differential expression due to ethanol

########### IF YOU NEED TO DOWNLOAD DESeq2 #############
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")#upgrade biocLite if you need to
# biocLite("DESeq")
# biocLite("DESeq2")
# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")
######################################################

library(DESeq2)
library(plotrix)

#----------------- UPLOAD DATA -----------------
#set working direcotry
source('deseq/zebrafish_RNAseq_functions.R')


#upload counts table
lnames=load('datasets/DEseq_inputs.Rdata') #file output from get_variance_stabilized_counts.R. Still includes outliers. Have not been filtered by basemean or mean readcount.
lnames


#look at the head of the file
head(counts)
dim(counts)


#get new totals
tots = apply(counts, 2, sum)
mean(tots)
std.error(tots)



#remove the 6-hour samples because they lack ethanol treatments
counts = counts[,-grep('NoE.6', colnames(counts))]
coldata=coldata[-grep('NoE.6', coldata$sample.names),]
dim(counts)
dim(coldata) #5 samples removed


#how many genes have average count higher than 3?
mns = apply(counts, 1, mean)
table(mns > 3)

#note categorical treatment of time
class(coldata$time)

#===============================================


#---------- BUILD A DESeq INPUT TABLES ----------
#one for the full model
#with factors for treatment, time, sequencing job, and treatment-time interaction
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~treatment + time + seqjob))
#==================================================


################ FULL MODEL ETHANOL EFFECT ###################
dds.eth=DESeq(ddsHTSeq) #this returns a single set of p values indicating the significance of the increase in likelihood gained by including the ethanol factor.
resultsNames(dds.eth)
res.eth=results(dds.eth, contrast = c('treatment', 'e', 'c'), independentFiltering=F)       #although a single set of p values is returned, the test will still have multiple possible contrasts which you can access, but these will only includence the log2 fold differences. Here we pick c('treatment', 'e', 'c') so that log2 fold differences describe ethanol treatment, with positive values indicated upregulation in ethanol treated samples and negative values downregulation in ethanol treated samples.
res.eth = res.eth[order(res.eth$pvalue),]
head(res.eth)
summary(res.eth)
volcano_plot(res.eth, MAIN='Full Dataset\nEthanol Effect')
res.names = merge_gene_names(data.frame(res.eth), sort.column='pvalue')
dim(res.names)
write.table(res.names, "results/ethanol_LRT_results.tsv", quote = F, sep = "\t")
save(res.eth, file='deseq/ethanol_full_LRT_results.Rdata')

#output significant genes
#these can be used for GO enrichment on the Gene ontology website
sig = na.omit(res.eth)
sig = sig[sig$padj < .01,]
head(sig)
nrow(sig)
write.table(sig, file = "results/significant_LRT_genes_ethanol.txt", quote = F, row.names = T)

#you can look at the actual difference in counts
for (g in rownames(sig[1:10,])){
	# pdf(file = paste(g, "counts_figure.pdf", sep = "_"))
	d <- plotCounts(dds.eth, gene=g, intgroup="treatment", returnData=F)
	# dev.off()
}

#for an individual gene of intersest
plotCounts(dds.eth, gene=g, intgroup="treatment", returnData=F)

###################################
####### output for GO-MWU #########
###################################
write.gomwu.input(res.eth, 'goMWU/ethanol_goMWU_input.csv')
###################################
###################################

#set working direcotry
save(coldata, counts, dds.eth, res.eth, file = 'deseq/deseq_objects.Rdata')
lnames = load('deseq/deseq_objects.Rdata')
lnames
####################################################

#--------- TEST FOR ETHANOL EFFECTS AT INDIVIDUAL TIME POINTS ---------------
#set up reduced coldata sets
edat = coldata[coldata$time==8,]
tdat = coldata[coldata$time==10,]
fdat = coldata[coldata$time==14,]


#### TEST WITHIN 8HR TIMEPOINT
#subset the counts
ecounts = counts[,colnames(counts) %in% edat$sample.names]
sum(colnames(ecounts) == edat$sample.names) == ncol(ecounts)

#set up DESeq input matrix
ddsHTSeq<-DESeqDataSetFromMatrix(ecounts, colData = edat, design = formula(~treatment+seqjob))

#run DESeq
dds.e=DESeq(ddsHTSeq)

#extract results
resultsNames(dds.e)
res.e =results(dds.e, contrast = c('treatment', 'e', 'c'), independentFiltering=F)
res.e = res.e[order(res.e$pvalue),]
head(res.e)
summary(res.e)
volcano_plot(res.e, MAIN='Ethanol Effect 8hr')

#write out
enames = merge_gene_names(data.frame(res.e), sort.column='pvalue')
save(dds.e, file='results/ethanol_8hr_LRT_results.Rdata')
write.table(enames, "results/ethanol_8hr_LRT_results.tsv", quote = F, sep = "\t")


#### TEST WITHIN 10HR TIMEPOINT
#subset the counts
tcounts = counts[,colnames(counts) %in% tdat$sample.names]
sum(colnames(tcounts) == tdat$sample.names) == ncol(tcounts)
#set up DESeq input matrix
ddsHTSeq<-DESeqDataSetFromMatrix(tcounts, colData = tdat, design = formula(~treatment+seqjob))
#run DESeq
dds.t=DESeq(ddsHTSeq)
#extract results
resultsNames(dds.t)
res.t =results(dds.t, contrast = c('treatment', 'e', 'c'), independentFiltering=F)
res.t = res.t[order(res.t$pvalue),]
head(res.t)
summary(res.t)
volcano_plot(res.t, MAIN='Ethanol Effect 10hr')

#write out
tnames = merge_gene_names(data.frame(res.t), sort.column='pvalue')
save(dds.t, file='results/ethanol_10hr_LRT_results.Rdata')
write.table(tnames, "results/ethanol_10hr_LRT_results.tsv", quote = F, sep = "\t")

#### TEST WITHIN 14HR TIMEPOINT
#subset the counts
fcounts = counts[,colnames(counts) %in% fdat$sample.names]
sum(colnames(fcounts) == fdat$sample.names) == ncol(fcounts)
#set up DESeq input matrix
ddsHTSeq<-DESeqDataSetFromMatrix(fcounts, colData = fdat, design = formula(~treatment))
#run DESeq
dds.f=DESeq(ddsHTSeq)
#extract results
resultsNames(dds.f)
res.f =results(dds.f, contrast = c('treatment', 'e', 'c'), independentFiltering=F)
res.f = res.f[order(res.f$pvalue),]
head(res.f)
summary(res.f)
volcano_plot(res.f, MAIN='Ethanol Effect 14hr')

#write out
fnames = merge_gene_names(data.frame(res.f), sort.column='pvalue')
save(dds.f, file='results/ethanol_14hr_LRT_results.Rdata')
write.table(fnames, "results/ethanol_14hr_LRT_results.tsv", quote = F, sep = "\t")

#--------- SAVE/LOAD ---------------
directory<-"~/gitreps/zebrafish_early_ethanol_RNASeq"
setwd(directory)
source('scripts/zebrafish_RNAseq_functions.R')
# save(dds.e, res.e, dds.t, res.t, dds.f, res.f, file="results/time_split_ethanol_tests.Rdata")
lnames = load('results/ethanol_full_LRT_results.Rdata');lnames
lnames = load('results/time_split_ethanol_tests.Rdata');lnames


#--------- HEATMAPS ---------------
library(pheatmap)
#load variance stabilized counts
lnames=load("~/gitreps/zebrafish_early_ethanol_RNASeq/datasets/raw_rld.Rdata")#variance stabilized counts for job2 only. Output from script get_variance_stabilized_counts.R
lnames = load("~/gitreps/zebrafish_early_ethanol_RNASeq/datasets/outliers.Rdata")
rld.df=rld.df[,!colnames(rld.df) %in% outlierNames]
head(rld.df)
dim(rld.df)

#subset for significant ethanol genes
head(res.eth)
CUT=0.05
sig=res.eth[!is.na(res.eth$padj) & res.eth$padj<CUT,]
dim(sig)
sig.rld=rld.df[rownames(rld.df) %in% rownames(sig),]

#get names
names = get_gene_names(rownames(sig.rld))
labs = c()
for (i in 1:nrow(names)){
	n=names$description[i]
	if (n ==''){
		n=names$external_gene_name[i]
	}
	labs=append(labs, strsplit(n, '[Source', fixed=T)[[1]][1])	
}
labs

#plot the heatmap
# heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=3)(100)
pheatmap(sig.rld,cluster_cols=T, cluster_rows=T,border_color=NA,clustering_distance_rows="correlation", labels_row=labs)



#--------- VOLCANO PLOTS ---------------#
source('scripts/zebrafish_RNAseq_functions.R')
XLIM=c(-4.2,4.2)
YLIM=c(0, 35)
volcano_plot(data.frame(res.eth), MAIN='Ethanol Effect\nAll Timepoints', XLIM=XLIM)
volcano_plot(res.e, MAIN='Ethanol Effect\n8hr', XLIM=XLIM)
volcano_plot(res.t, MAIN='Ethanol Effect\n10hr', XLIM=XLIM)
volcano_plot(res.f, MAIN='Ethanol Effect\n14hr', XLIM=XLIM)

#ggplot
NAME=F
top=10

a=ggvolcano_plot(res.eth, XLIM=XLIM, YLIM=F, addNames = NAME, topN = top, MAIN='Ethanol Effect', submain='All Timepoints', xshift=0.4, yshift=0.8)
e=ggvolcano_plot(res.e, XLIM=XLIM, YLIM=F, addNames = NAME, topN = top, MAIN='Ethanol Effect', submain='8 hour', xshift=c(0.4), yshift=0.6)
t=ggvolcano_plot(res.t, XLIM=XLIM, YLIM=F, addNames = NAME, topN = top, MAIN='Ethanol Effect', submain='10 hour', xshift=c(-0.4, 0.4), yshift=0.6)
f=ggvolcano_plot(res.f, XLIM=XLIM, YLIM=F, addNames = NAME, topN = top, MAIN='Ethanol Effect', submain='14 hour', xshift=0.4, yshift=0.4)
write.table(a, file='~/Desktop/topAll.tsv')
write.table(e, file='~/Desktop/top8.tsv')
write.table(t, file='~/Desktop/top10.tsv')
write.table(f, file='~/Desktop/top14.tsv')


#plot with color coding for overlapping significant between datasets
s1 = rownames(res.eth)[res.eth$padj<0.1 & !is.na(res.eth$padj)]
s2 = rownames(res.e)[res.e$padj<0.1 & !is.na(res.e$padj)]
s3 = rownames(res.t)[res.t$padj<0.1 & !is.na(res.t$padj)]
s4 = rownames(res.f)[res.f$padj<0.1 & !is.na(res.f$padj)]
s = s1[s1 %in% s2 & s1 %in% s3 & s1 %in% s4]
length(s)
sdf = data.frame(s)
rownames(sdf) = sdf$s
snames = merge_gene_names(sdf)
snames$s<-NULL
write.table(snames, "~/junk/sigAllTimepoints.txt", row.names=F, quote=F, sep="\t")



res.eth$inc = rownames(res.eth) %in% s
res.eth$sig = res.eth$padj < 0.1
sub.eth = data.frame(res.eth[rownames(res.eth) %in% s,])
g = ggplot(data= data.frame(res.eth), aes(x=log2FoldChange, y=-log10(pvalue), colour=sig)) + 
	scale_colour_manual(values=c('black', 'red')) +
	geom_point(alpha=0.4, size=1.75) + 
	
	geom_point(data=sub.eth, aes(x=log2FoldChange, y=-log10(pvalue)), color='blue') +
	xlab("log2 fold difference") + 
	ylab("-log10 p-value") +
	theme_bw() 
plot(g)




require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
sig.col='red'
ns.col='black'
addNames = TRUE
topN = 10
deseq.res=na.omit(data.frame(res.eth))
deseq.res$threshold = deseq.res$padj < 0.1
tolab=deseq.res[1:topN,]
sig = merge_gene_names(tolab)

##Construct the plot object
g = ggplot(data= deseq.res, aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) + scale_colour_manual(values=c(ns.col, sig.col)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(XLIM) + ylim(YLIM) +
  xlab("log2 fold change") + ylab("-log10 p-value") + theme_bw() + guides(colour=FALSE) + ggtitle('Ethanol Effect', subtitle='All Timepoints')

if (addNames){
	g + geom_text(data=sig, aes(x=sig$log2FoldChange, y=-log10(sig$pvalue),
	                     label=sig$external_gene_name, size=1.2), colour="black") + guides(size=FALSE)
}







#--------- COMPARE ETHANOL EFFECTS ACCROSS TIME POINTS ---------------
source('scripts/venn_diagram_functions.R')
CUT=0.1
get.sig=function(res, TIME){
	sig=data.frame(res[!is.na(res$padj) & res$padj<=CUT,])
	sig$time=TIME
	x=nrow(sig)
	print(paste(x, 'significant genes'))
	sig=sig[order(sig$pvalue),]
	return(sig)
}


#subset for significant genes
sig.all=get.sig(res.eth, TIME='all')
sig.e=get.sig(res.e, TIME=8)
sig.t=get.sig(res.t, TIME=10)
sig.f=get.sig(res.f, TIME=14)


#write these out for easy access
write.table(merge_gene_names(sig.all), file='results/sig_eth_full_dataset.tsv', quote=F, sep='\t')
write.table(merge_gene_names(sig.e), file='results/sig_eth_8hr.tsv', quote=F, sep='\t')
write.table(merge_gene_names(sig.t), file='results/sig_eth_10hr.tsv', quote=F, sep='\t')
write.table(merge_gene_names(sig.f), file='results/sig_eth_14hr.tsv', quote=F, sep='\t')




#LOOK AT PAIRWISE VENN DIAGRAMS
pairwiseVenn(sig.all, sig.e, 'Full Dataset', '8hr')
pairwiseVenn(sig.all, sig.t, 'Full Dataset', '10hr')
pairwiseVenn(sig.all, sig.f, 'Full Dataset', '14hr')
pairwiseVenn(sig.e, sig.t, '8hr', '10hr')
pairwiseVenn(sig.e, sig.f, '8hr', '14hr')
pairwiseVenn(sig.t, sig.f, '10hr', '14hr')
pairwiseVenn(sig.t, sig.f, '10hr', '14hr')

#LOOK AT 4-WAY VENN DIAGRAM
s1=rownames(sig.all)
s2=rownames(sig.e)
s3=rownames(sig.t)
s4=rownames(sig.f)
CATEGORY=c('Full', '8hr', '10hr', '14hr')
intersection.list = fourway_venn(s1, s2, s3, s4, CATEGORY=rep('',4))
print(intersection.list)


#output the intersections
for (i in 1:length(intersection.list)){
	name=names(intersection.list)[i]
	x=get_gene_names(intersection.list[i])
	filename=paste(name, '_intersect_genes.tsv', sep='')
	print(filename)
	print(x)
	write.table(x, file=paste('results/', filename, sep=''))
}


#output set of genes significant in any ethanol test
any.eth = unique(c(s1, s2, s3, s4))
write.table(any.eth, file='results/sig_eth_any_test.tsv', quote=F, row.names=F, sep="\t")



#LOOK AT GENES IN INTERESTING CATEGORIES
source('scripts/zebrafish_RNAseq_functions.R')
sa = merge_gene_names(sig.all, sort.column='pvalue')
se = merge_gene_names(sig.e)
st = merge_gene_names(sig.t)
sf = merge_gene_names(sig.f)


