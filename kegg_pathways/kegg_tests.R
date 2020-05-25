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



##!!! Note! Tried running this with same.dir=TRUE, and got lots of keggs that were 'lower'
#but when I looked at the actual expression values, they did not seem convincing at all
#keeping the gage stuff here for reference, but dont' really trust it.



# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
rm(list=ls())
library(dplyr)
library(DESeq2)

# Import DESeq2 results to use
lnames = load('deseq/ethanol_full_LRT_results.Rdata')  #gerated by DEseq_early_ethanol.R 
lnames
res=res.eth


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
library(gage)


#the package gageData does not have the zebrafish Keggs, so they need to be done manually
#upload zebrafish keggs retrieved from here: http://rest.kegg.jp/link/dre/pathway
keggs = read.table("metadata/zebrafish_keggs.tsv", header = T, colClasses = c("character", "character")) #get this with KeggAPI (http://rest.kegg.jp/link/dre/pathway)
keg.names = read.table("metadata/zebrafish_kegg_names.tsv", header = T, sep = "\t", colClasses=c("character", "character")) #get this with keggAPI (http://rest.kegg.jp/list/pathway/dre)
head(keggs)
head(keg.names)
#build the kegg list
name.set = keg.names$keg
kegg.sets.dr = lapply(name.set, function(x) keggs$gene[keggs$keg == x])
names(kegg.sets.dr) = name.set
head(kegg.sets.dr, 3)
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)



###For experimentally derived gene sets, GO term groups, etc, coregulation is commonly the case, hence same.dir = TRUE (default); In KEGG, BioCarta pathways, genes frequently are not coregulated, hence it could be informative to let same.dir = FALSE. Although same.dir = TRUE could also be interesting for pathways.
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.dr, same.dir=FALSE) #note the same.dir argument here to decide the 'signedness' of the test

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
# p.type = 'p.val'
p.type = 'padj'
CUT = 0.05


head(keggres)
#get dataframes for up and downregulated keggs (greater and less)
gage_res = na.omit(data.frame(keggres$greater))
gage_sig = gage_res[gage_res[,p.type] <= CUT,]
nrow(gage_sig)


#save objects for plotting selected keggs
save(foldchanges, res, keggs, file='kegg_pathways/pathview_input.Rdata')


# test for enrichment with fisher's exact tests ---------------------------
kdat = keggs %>% 
  dplyr::rename(entrez = gene) %>% 
  left_join(data.frame(res), by = 'entrez') %>% 
  as_tibble()

fdr_alpha = 0.05

fisher_results = data.frame()
for (k in unique(kdat$kegg)){
  ksub = kdat %>% 
    mutate(in_kegg = kegg == k,
           sig_up = log2FoldChange > 0 & padj < fdr_alpha,
           sig_down = log2FoldChange < 0 & padj < fdr_alpha)
  
  fisher_up = fisher.test(x = ksub$in_kegg,
                          y = factor(ksub$sig_up),
                          alternative = 'greater')
  fisher_down = fisher.test(x = ksub$in_kegg,
                            y = factor(ksub$sig_down),
                            alternative = 'greater')
  up_res = c(k, 'upregulation', fisher_up$estimate, fisher_up$p.value)
  down_res = c(k, 'downregulation', fisher_down$estimate, fisher_down$p.value)
  fisher_results = rbind(fisher_results, rbind(up_res, down_res))
}

#format the results
colnames(fisher_results) = c('kegg', 'direction', 'odds_ratio', 'pvalue')
fisher_results$pvalue = as.numeric(as.character(fisher_results$pvalue))
fisher_results$padj = p.adjust(fisher_results$pvalue, method='BH')
fisher_results = fisher_results %>% 
  as_tibble()

#check significant
fisher_sig = fisher_results %>% 
  filter(padj < 0.05)
fisher_sig
keggresids = fisher_sig$kegg

# Define plotting function for applying later
#set the directory where you want to output the results
# setwd("gage/ethanol_results")
LOW='dodgerblue'
MID='grey'
HIGH='red'
scale.limit=1
NODE.SUM = 'max.abs'
sig_entrez = data.frame(res) %>% 
  filter(padj < 0.05) %>% 
  pull(entrez)
sig_fold_changes = foldchanges[sig_entrez]

# plot multiple pathways (plots saved to disk and returns a throwaway list object)

#plot pdfs
setwd('./kegg_pathways/significant_pathview_figures/')
tmp = sapply(keggresids,
             function(pid) 
               pathview(gene.data=sig_fold_changes,
                        pathway.id=pid,
                        species="dre",
                        node.sum= NODE.SUM,
                        low= LOW,
                        mid=MID,
                        high=HIGH,
                        limit=list(gene=scale.limit,
                                   cpd=scale.limit),
                        kegg.native=FALSE,
                        same.layer=F)
)

#plot pngs
tmp = sapply(keggresids, function(pid) pathview(gene.data=sig_fold_changes, pathway.id=pid, species="dre", node.sum= NODE.SUM, low= LOW, mid=MID, high=HIGH, limit=list(gene=scale.limit,cpd=scale.limit), kegg.native=TRUE))
setwd('../../')
fisher_sig

# assemble keggs and genes for checking -----------------------------------

#assemble a single dataframe
kegg_vector =  unlist(kegg.sets.dr[keggresids])
gene_df = keggs %>% 
  dplyr::rename(entrez = gene) %>% 
  filter(kegg %in% keggresids) %>% 
  dplyr::inner_join(data.frame(res), by = 'entrez') %>% 
  as_tibble()
unique(gene_df$kegg)

# write out dataframes for each selected kegg -----------------------------

for (i in 1:length(keggresids)){
  my_kegg = keggresids[i]
  print(paste(my_kegg, '...', sep=''))
  out_name = paste(my_kegg, 'gene_results.csv', sep='_')
  out_path = paste('kegg_pathways/significant_pathview_figures/', out_name, sep='')
  kdat = gene_df %>% 
    dplyr::filter(kegg==my_kegg)
  kdat %>% 
    write_csv(out_path)
}

