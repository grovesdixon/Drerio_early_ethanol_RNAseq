#gage_ethanol.R
#Groves Dixon
#10/26/16
#based on tutorial 'RNA-seq differential expression & pathway analysis with Sailfish, DESeq2, GAGE, and Pathview'
#found here: https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/

########### IF YOU NEED TO DOWNLOAD GAGE #############
# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")#upgrade biocLite if you need to
# biocLite("gage")
# biocLite(c("pathview", "gage", "gageData", "GenomicAlignments","TxDb.Hsapiens.UCSC.hg19.knownGene"))
# ######################################################

# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)
setwd('/Users/grovesdixon/gitreps/Drerio_early_ethanol_RNAseq/')

# Import DESeq2 results to use
lnames = load('deseq/ethanol_full_LRT_results.Rdata')  #gerated by DEseq_early_ethanol.R 
lnames
res=res.eth



#----- OPTIONALLY SWITCH TO SIGNIFICANT ONLY! -----#
#Choose if you want to use significant genes only
#If TRUE, then only significant genes will be color coded in pathway figures
#If FALSE, then all genes with expression data will be color coded
sigOnly=FALSE

before = nrow(res)
if (sigOnly){
  setwd('./kegg_pathways/sig_only_results')
  print('RUNNGING FOR ONLY SIGNIFICANT GENES')
  res=res[!is.na(res$padj) & res$padj < 0.1,]
  after=nrow(res)
  print(paste(c('Building plots for', after, 'significant genes'), collapse=' '))
} else{
  setwd('./kegg_pathways/ethanol_results')
  print('Building pathways for all expressed genes')
}
#--------------------------------------------------#

#sort and summarize
res = res[order(res$pvalue),]
summary(res)


###"Since we mapped and counted against the Ensembl annotation, our results only have information about Ensembl gene IDs. But, our pathway analysis downstream will use KEGG pathways, and genes in KEGG pathways are annotated with Entrez gene IDs. I wrote an R package for doing this offline the dplyr way (https://github.com/stephenturner/annotables), but the canonical Bioconductor way to do it is with the AnnotationDbi and organism annotation packages. Here we’re using the organism package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”) using Entrez Gene IDs (“eg”) as primary keys. To see what all the keys are, use the columns function."--tutorial
#upload libraries for getting annotations for ensemble ids
library("AnnotationDbi")
# library("org.Hs.eg.db") #use this one for human
library("org.Dr.eg.db") #use this one for Zebrafish
columns(org.Dr.eg.db)

#grab annotations and add them to the results dataframe
res$symbol = mapIds(org.Dr.eg.db,
                     keys=row.names(res), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez = mapIds(org.Dr.eg.db,
                     keys=row.names(res), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first") #here we get the Entrez ID, which is used to link genes to keggs
res$name =   mapIds(org.Dr.eg.db,
                     keys=row.names(res), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

head(res, 10)



#now use gage on the data
library(pathview)
library(gage)


#the package gageData does not have the zebrafish Keggs, so they need to be done manually
#upload zebrafish keggs retrieved from here: http://rest.kegg.jp/link/dre/pathway
keggs = read.table("datasets/zebrafish_keggs.tsv", header = T, colClasses = c("character", "character")) #get this with KeggAPI (http://rest.kegg.jp/link/dre/pathway)
keg.names = read.table("datasets/zebrafish_kegg_names.tsv", header = T, sep = "\t", colClasses=c("character", "character")) #get this with keggAPI (http://rest.kegg.jp/list/pathway/dre)
head(keggs)
head(keg.names)
#build the kegg list
name.set = keg.names$keg
kegg.sets.dr = lapply(name.set, function(x) keggs$gene[keggs$keg == x])
names(kegg.sets.dr) = name.set
head(kegg.sets.dr, 3)


### The gage() function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs.
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)



###For experimentally derived gene sets, GO term groups, etc, coregulation is commonly the case, hence same.dir = TRUE (default); In KEGG, BioCarta pathways, genes frequently are not coregulated, hence it could be informative to let same.dir = FALSE. Although same.dir = TRUE could also be interesting for pathways.
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.dr, same.dir=F) #note the same.dir argument here to decide the 'signedness' of the test

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

##Now, let’s process the results to pull out the top 5 upregulated pathways, then further process that just to get the IDs. We’ll use these KEGG pathway IDs downstream for plotting.
# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways


# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids


#perform FDR and subset for significant
#select whether you want to use adjusted or normal p values and the cutoff from printing out path figures
p.type = 'p.val'
p.type = 'padj'
CUT = 0.05


head(keggres)
#get dataframes for up and downregulated keggs (greater and less)
greater = na.omit(data.frame(keggres$greater))
less = na.omit(data.frame(keggres$less))
total.tests = nrow(greater) + nrow(less)
greater$padj = p.adjust(greater$p.val, n=total.tests, method='BH')
less$padj=p.adjust(less$p.val, n=total.tests, method='BH')
all=rbind(greater, less)
#subset for sig paths
sig.greater = greater[greater[,p.type] <= CUT,]
sig.less = less[less[,p.type] <= CUT,]
sig.greater
sig.less
sig.all = rbind(sig.greater, sig.less)

keggrespathways =append(rownames(sig.greater), rownames(sig.less))
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids



#detach dply
detach("package:dplyr", unload=TRUE)

# # Define plotting function for applying later
# #set the directory where you want to output the results
# plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="dre", new.signature=FALSE)
# 
# 
# 
# # plot multiple pathways (plots saved to disk and returns a throwaway list object)
# #default settings for testing:
# tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="dre"))
# write.table(sig.all, file="kegg_stats.tsv", sep = "\t", quote = F) #write out the results
# 
# #plot pngs
# #modified
# LOW='dodgerblue'
# MID='grey'
# HIGH='red'
# scale.limit=1
# tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="dre", node.sum='sum', low= LOW, mid=MID, high=HIGH, limit=list(gene=scale.limit,cpd=scale.limit), kegg.native=TRUE))
# 
# #OUTPUT A SELECTED SET OF KEGG PATHWAYS
# Ribosome  = dre03010
# Apoptosis = dre04210
# Wnt       = dre04310 
# Hedgehod  = dre04340
# PI3K/AKT  = dre04151

# mTOR signaling pathway (dre04150)
# oxidative phosphorylation (dre00190)
# DNA replication (dre03030)
# cell cycle (dre04110)
# p53 signaling (dre04115)
# notch signaling (dre04330) 
# TGF-beta signaling (dre04340)

#enter as vectors
keggresids=c('dre03010',
             'dre04210',
             'dre04310',
             'dre04340',
             'dre04150',
             'dre00190',
             'dre03030',
             'dre04110',
             'dre04115',
             'dre04330',
             'dre04340') #alfire's choices
knames = c("Ribosome",
           "Apoptosis",
           "Wnt",
           "Hedgehod",
           "mTOR",
           "oxidative phosphorylation",
           "DNA replication",
           "cell cycl",
           "p53",
           "notch",
           "TGF-beta")


#or load them individually
keggresids = c("dre04340")  #Hedgehog
knames = c('Hedgehod')      #Hedgehog

keggresids = c("dre04310")  #wntPCP
knames = c('Wnt')      #wntPCP

#run pathview as above
dev.res = all[rownames(all) %in% keggresids,]
dev.res

# Define plotting function for applying later
#set the directory where you want to output the results
# setwd("gage/ethanol_results")
LOW='dodgerblue'
MID='grey'
HIGH='red'
scale.limit=1
NODE.SUM = 'max.abs'

# plot multiple pathways (plots saved to disk and returns a throwaway list object)

#plot pdfs
tmp = sapply(keggresids,
             function(pid) 
               pathview(gene.data=foldchanges,
                        pathway.id=pid,
                        species="dre",
                        node.sum= NODE.SUM,
                        low= LOW,
                        mid=MID,
                        high=HIGH,
                        limit=list(gene=scale.limit,cpd=scale.limit),
                        kegg.native=FALSE,
                        same.layer=F)
             )

#plot pngs
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="dre", node.sum= NODE.SUM, low= LOW, mid=MID, high=HIGH, limit=list(gene=scale.limit,cpd=scale.limit), kegg.native=TRUE))


#view genes associated with the KEGGS
head(kegg.sets.dr)
keggGenes = unlist(kegg.sets.dr[keggresids])
res.df = data.frame(res)
keggResults = res.df[res.df$entrez %in% keggGenes,]

#look at all instances for an individual gene
target.gene = 'shhb' #note these match with the pdf pathview
target.gene = 'shha' #note these match with the pdf pathview
targetResults = keggResults[keggResults$symbol==target.gene,]
targetResults


library(tidyverse)
keggResults %>% 
  filter(grepl('rac', symbol))

keggResults[grepl('rac', keggResults$symbol),] %>% 
  write_tsv(path='~/junk/racs.tsv')
