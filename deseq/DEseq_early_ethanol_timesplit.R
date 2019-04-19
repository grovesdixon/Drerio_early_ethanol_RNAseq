#DEseq_early_ethanol_timesplit.R

lnames = load('deseq/deseq_objects.Rdata')
lnames

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
save(dds.e, res.e, file='deseq/ethanol_8hr_LRT_results.Rdata')
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
save(dds.t,res.t file='deseq/ethanol_10hr_LRT_results.Rdata')
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
save(dds.f, res.f, file='deseq/ethanol_14hr_LRT_results.Rdata')
write.table(fnames, "results/ethanol_14hr_LRT_results.tsv", quote = F, sep = "\t")


#--------- VOLCANO PLOTS ---------------#
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









#--------- COMPARE ETHANOL EFFECTS ACCROSS TIME POINTS ---------------
source('deseq/venn_diagram_functions.R')
CUT=0.1

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

#LOOK AT 3-WAY VENN
s1=rownames(sig.all)
s2=rownames(sig.e)
s3=rownames(sig.t)
s4=rownames(sig.f)
intersection.list = threeway_venn(s2, s3, s4, CATEGORY=c('8hr', '10hr', '14hr'))
print(intersection.list)


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
sa = merge_gene_names(sig.all, sort.column='pvalue')
se = merge_gene_names(sig.e)
st = merge_gene_names(sig.t)
sf = merge_gene_names(sig.f)


