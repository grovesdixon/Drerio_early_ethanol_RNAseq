#figure_plotting.R
#plot the figures for publication

library(tidyverse)
library(dplyr)
library(cowplot)


setwd('~/gitreps/Drerio_early_ethanol_RNAseq/')
lnames=load("datasets/large_ignored/initialize_countsImage1.Rdata")
source("deseq/zebrafish_RNAseq_functions.R")

# main heatmap ------------------------------------------------------------


x=colnames(rld.df)
y=sub("NoE", "C", x)
z=sub("h.", "", y)
l = substr(z, 1,nchar(z)-1)
l2=sub(".", "", l, fixed=T)
pheatmap(cor(rld.df, method = 'spearman'), labels_row=l2, labels_col=l2, treeheight_row=0, treeheight_col=0, number_color='blue')



# pcas --------------------------------------------------------------------

ptSIZE=3

#developmental stage
devPca = rld.df %>% 
  mod.plotPCA.df(coldat = coldata, intgroup = 'time', main = 'Timepoint', legendTitle='Hpf', SIZE = ptSIZE)
# devPca

#batch pc1 and pc2
batchPca = rld.df %>% 
  mod.plotPCA.df(coldat = coldata,
                 intgroup = 'seqjob',
                 main = "Batch",
                 legendTitle='Batch',
                 SIZE = ptSIZE
  )
# batchPca 


#ethanol pc1 and pc2
ethPca = rld.df %>% 
  mod.plotPCA.df(coldat = coldata %>% 
                   mutate(t2=if_else(treatment=='c',
                                     'Control',
                                     'Ethanol')),
                 intgroup = 't2',
                 main = "Ethanol Treatment",
                 legendTitle='Treatment', 
                 SIZE = ptSIZE
  )
# ethPca

#ethanol pc1 and pc2
ethPca89 = rld.df %>% 
  mod.plotPCA.df(pc1 = 8,
                 pc2 = 9,
                 coldat = coldata %>% 
                   mutate(t2=if_else(treatment=='c',
                                     'Control',
                                     'Ethanol')),
                 intgroup = 't2',
                 main = "Ethanol Treatment",
                 legendTitle='Treatment', 
                 SIZE = ptSIZE+2
  )
# ethPca89

plot_grid(devPca, batchPca, ethPca, nrow=3)
quartz()
plot(ethPca89)

# volcano plots -----------------------------------------------------------

#individual
NAME=F
top=10
source('deseq/zebrafish_RNAseq_functions.R')
ll=load('deseq/deseq_objects.Rdata')
ll=load('deseq/ethanol_8hr_LRT_results.Rdata')
ll=load('deseq/ethanol_10hr_LRT_results.Rdata')
ll=load('deseq/ethanol_14hr_LRT_results.Rdata')
XLIM=c(-4.2,4.2)
YLIM=c(0, 35)
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
write.table(snames, "results/sigAllTimepoints.txt", row.names=F, quote=F, sep="\t")
res.eth$inc = rownames(res.eth) %in% s
res.eth$sig = res.eth$padj < 0.1
sub.eth = data.frame(res.eth[rownames(res.eth) %in% s,])
g = data.frame(res.eth) %>% 
  mutate(sig=if_else(is.na(sig),
                     FALSE,
                     sig)) %>% 
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue), colour=sig)) + 
  scale_colour_manual(values=c('black', 'red')) +
  geom_point(alpha=0.4, size=1.75) + 
  
  geom_point(data=sub.eth, aes(x=log2FoldChange, y=-log10(pvalue)), color='blue') +
  xlab("log2 fold difference") + 
  ylab("-log10 p-value") +
  guides(color=guide_legend(title="FDR < 0.1"))
plot(g)



# venn diagram ------------------------------------------------------------
source('deseq/venn_diagram_functions.R')
CUT=0.1
sig.e=get.sig(res.e, TIME=8)
sig.t=get.sig(res.t, TIME=10)
sig.f=get.sig(res.f, TIME=14)
s2=rownames(sig.e)
s3=rownames(sig.t)
s4=rownames(sig.f)
intersection.list = threeway_venn(s2, s3, s4, CATEGORY=c('8hr', '10hr', '14hr'))
print(intersection.list)



# plot time modules ----------------------------------------------------------------
ll=load('results/timeMods.Rdata')
ll=load('datasets/large_ignored/raw_rld.Rdata')
rld.df=data.frame(rld.df)
modules = timeMods$module
geneModuleMembership$gene = rownames(geneModuleMembership)
head(geneModuleMembership)
hubGenes = c()
for (m in modules){
  sub=geneModuleMembership[order(abs(geneModuleMembership[,m])),c(m,'gene')]
  hub=sub[nrow(sub),'gene']
  hubGenes = append(hubGenes, hub)
}
hdat = data.frame('module'=modules,
              'hub' = hubGenes,
              row.names=hubGenes, stringsAsFactors=F) %>% 
  merge_gene_names()
hdat

#plot each of the time module hub genes for ethanol and controls
plotList=list()
for (i in 1:nrow(hdat)){
  print(paste(i,'...',sep=''))
  row=hdat[i,]
  GENE=row[,'hub']
  NAME = row[,'external_gene_name']
  MODULE = sub('ME', '', row[,'module'])
  if(is.na(NAME)){
    NAME='none'
  }
  tdat = rld.df %>% 
    mutate(gene=rownames(rld.df)) %>% 
    filter(gene==GENE) %>% 
    dplyr::select(-gene) %>% 
    t() %>% 
    data.frame()
  colnames(tdat) = 'ge'
  s=sub('NoE.', 'C', rownames(tdat), fixed=T)
  s=sub('E.', 'E', s, fixed=T)
  s=sapply(s, function(x) strsplit(x, '.', fixed=T)[[1]][1])
  tdat$treat = substr(s, 1,1)
  tdat$time=sub('h', '', substr(s, 2,5)) %>% as.numeric()
  plt=tdat %>% 
    group_by(treat, time) %>% 
    ggplot(aes(x=time, y=ge, color=treat, fill=treat)) +
    geom_jitter(width=0.2) +
    geom_smooth(se=T) +
    labs(x='hpf', y='Expression level',
         subtitle=paste(GENE,NAME,sep='\n'),
         title=MODULE) +
    theme(plot.title = element_text(colour = 'black'),
          legend.title = element_blank()) +
    scale_x_continuous(breaks=c(6,8,10,14))
  plotList[[i]]=plt
}

plot_grid(plotlist = plotList, nrow=2)



# plot gpc4 ---------------------------------------------------------------

#pick the gene
g = c("ENSDARG00000015472")

#pick the time point
lnames = load('deseq/ethanol_full_LRT_results.Rdata') #All 
resultsNames(dds.eth)
myPlotCount(dds.eth, g, 'All timepoints')


#run for each individual timepoint
ll = load('deseq/ethanol_8hr_LRT_results.Rdata');print(ll)
ll = load('deseq/ethanol_10hr_LRT_results.Rdata');print(ll)
ll = load('deseq/ethanol_14hr_LRT_results.Rdata');print(ll)
full = myPlotCount(dds.eth, g, 'All timepoints')
eight = myPlotCount(dds.e, g, '8hr') + theme(axis.title.y=element_blank())
ten = myPlotCount(dds.t, g, '10hr') + theme(axis.title.y=element_blank())
fort = myPlotCount(dds.f, g, '14 hr') + theme(axis.title.y=element_blank())
plot_grid(full, eight, ten, fort, nrow=1, label_x = 'Treatment', label_y='Normalized count', align='v', rel_widths=c(1,1,1,1))

