#figure_plotting.R
#plot the figures for publication

library(tidyverse)
library(cowplot)



lnames=load("datasets/large_ignored/initialize_countsImage1.Rdata")
source("./initialize_counts/multivariate_functions.R")


# main heatmap ------------------------------------------------------------


x=colnames(rld.df)
y=sub("NoE", "C", x)
z=sub("h.", "", y)
l = substr(z, 1,nchar(z)-1)
l2=sub(".", "", l, fixed=T)
pheatmap(cor(rld.df, method = 'spearman'), labels_row=l2, labels_col=l2, treeheight_row=0, treeheight_col=0, number_color='blue')



# pcas --------------------------------------------------------------------

#developmental stage
devPca = rld.df %>% 
  mod.plotPCA.df(coldat = coldata, intgroup = 'time', main = 'Time (hr)', legendTitle='Hpf')
devPca

#ethanol pc1 and pc2
ethPca = rld.df %>% 
  mod.plotPCA.df(coldat = coldata %>% 
                   mutate(t2=if_else(treatment=='c',
                                     'Control',
                                     'Ethanol')),
                 intgroup = 't2',
                 main = "Ethanol Treatment",
                 legendTitle='Treatment'
                 )
ethPca

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
                 legendTitle='Treatment'
  )
ethPca89


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
