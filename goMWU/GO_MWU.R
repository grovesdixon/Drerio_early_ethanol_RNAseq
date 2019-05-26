
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("./goMWU")


#SELECT WHICH DATASET TO INPUT
input = "ethanol_goMWU_input.csv" ;ABSVALUE=-log(0.05,10)    #deseq results for full dataset
input = "ethanol_up_mod_genes.csv";ABSVALUE=c(1e-5)   #upreulgated module from WGCNA
input = "ethanol_down_mod_genes.csv";ABSVALUE=c(1e-5) #downregulated module from WGCNA


# Edit these to match your data file names: 
goAnnotations="zebrafish_embl_go.tsv" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
	largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest=5,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results
quartz()
gomwuPlot(input,goAnnotations,goDivision,
	absValue=ABSVALUE,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
	level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
	level2=0.05, # FDR cutoff to print in regular (not italic) font.
	level3=0.01, # FDR cutoff to print in large bold font.
	txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
	treeHeight=0.5, # height of the hierarchical clustering tree
#	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.05,level2=0.01,level3=0.001.  




#------------ UPLOAD THE RESULTS FILE WITH STATS FOR EACH GO TERM ------------#
resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
go.res=read.table(resName, header = T)
go.res=go.res[order(go.res$pval),]
head(go.res, n=20)

#subset for significant GO terms
CUT=0.1
sigGo=go.res[round(go.res$p.adj, digits=1)<=CUT,]
head(sigGo)
nrow(sigGo)


#------------ LOOK AT GENES WITHIN GO TERMS ------------#
#upload the results for each genes output by the functions above
resName=paste(goDivision, input, sep='_')
gene.res=read.table(resName, header=T)
head(gene.res)

#select the GO term from your results file (sigGo above)
go="GO:0004984" #olfactory receptor activity
go='GO:0004930' #G-protein coupled receptor activity
go='GO:0007409;GO:0048812' #neuron projection morphogenesis

#subset for that GO term
go.genes = gene.res[gene.res$term == go, 'seq']
length(go.genes)

#get gene names
geneSet = as.character(go.genes)

#gather the names from ensemble using biomarRt
library("biomaRt")
embl = useMart("ensembl", dataset="drerio_gene_ensembl")
names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description'), filters = c('ensembl_gene_id'), values = geneSet, mart = embl)
head(names)

#named genes are here:
names


#------------ BUILD HEATMAP FOR A GIVEN GO TERM ------------#
library(pheatmap)


#you'll need the variance stabilized counts
lnames=load("~/git_repositories/zebrafish_early_ethanol_RNASeq/wgcna/wgcna-01_output_job2.Rdata")#variance stabilized counts for job2 only. Output from script get_variance_stabilized_counts.R
d=t(datExpr)
head(d)
dim(d)

#you'll also want the deseq results
lnames=load("~/git_repositories/zebrafish_early_ethanol_RNASeq/results/time_categorical_DESeq.Rdata")
x=data.frame(res.time)
s=x[!is.na(x$padj) & x$padj<0.0001,]
tail(s)
dim(s)
geneSet=as.character(rownames(s))

#subset the variance stabilized counts to include only the GO of interest
g=d[rownames(d) %in% geneSet,]
head(g)
dim(g)

#gather gene names from names above
# names$description[]

# colnames(labs)=c('ensembl_gene_id')
# labs=merge(labs, names, by = 'ensembl_gene_id', all.x=T)
# labs$ensembl_gene_id == rownames(g)

#plot the heatmap
pheatmap(g,cluster_cols=T,border_color=NA,clustering_distance_rows="correlation")









