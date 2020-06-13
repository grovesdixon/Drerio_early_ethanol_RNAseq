#figure_plotting.R
#plot the figures for publication
rm(list=ls())
lnames=load("datasets/raw_rld.Rdata")
source("deseq/zebrafish_RNAseq_functions.R")

# main heatmap ------------------------------------------------------------

library(pheatmap)
x=colnames(rld.df)
y=sub("NoE", "C", x)
z=sub("h.", "", y)
l = substr(z, 1,nchar(z)-1)
l2=sub(".", "", l, fixed=T)
pheatmap(cor(rld.df, method = 'spearman'), labels_row=l2, labels_col=l2, treeheight_row=0, treeheight_col=0, number_color='blue')



# pcas --------------------------------------------------------------------

ptSIZE=3

#developmental stage
coldata$time = factor(coldata$time, levels = c('6', '8', '10', '14'))


pca_df = build_pca(rld.df,
                   coldata,
                   ntop = 25000,
                   pcs = 9)

#make changes (tidyverse kills attr)
pca_df$PC1 = pca_df$PC1*-1
pca_df$etoh = if_else(pca_df$treatment == 'c',
                      'control',
                      'ethanol')

devPca = plot_rld_pca (pca_df,
              group_col = 'time',
              pc1 = 1,
              pc2 = 2,
              subtitle = "Age (hpf)",
              size = ptSIZE,
              legend_title=NULL,
              x_invert=1,
              legend_position = 'right',
              fix_coords = FALSE)


batchPca = plot_rld_pca (pca_df,
                       group_col = 'seqjob',
                       pc1 = 1,
                       pc2 = 2,
                       subtitle = "Batch",
                       size = ptSIZE,
                       legend_title=NULL,
                       x_invert=1,
                       legend_position = 'right',
                       fix_coords = FALSE)

ethPca = plot_rld_pca (pca_df,
                         group_col = 'etoh',
                         pc1 = 1,
                         pc2 = 2,
                         subtitle = "Treatment",
                         size = ptSIZE,
                         legend_title=NULL,
                         x_invert=1,
                         legend_position = 'right',
                         fix_coords = FALSE)

ethPca89 = plot_rld_pca (pca_df,
                       group_col = 'etoh',
                       pc1 = 8,
                       pc2 = 9,
                       subtitle = "Treatment",
                       size = ptSIZE,
                       legend_title=NULL,
                       x_invert=1,
                       legend_position = 'right',
                       fix_coords = FALSE)

plot_grid(devPca, batchPca, ethPca, ethPca89, nrow=2)


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

plot_grid(a, e, t, f, nrow=2)

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
  lims(x=c(-2,2)) +
  
  geom_point(data=sub.eth, aes(x=log2FoldChange, y=-log10(pvalue)), color='blue') +
  xlab(bquote(log[2]~"fold difference")) + 
  ylab('-'*log[10]~'p-value') +
  guides(color=guide_legend(title="FDR < 0.1")) +
  theme(legend.position = 'none')
plot(g)



# venn diagram ------------------------------------------------------------
#build a venn diagram of the overlapping significant genes from timepoints
library(limma)
library(ggforce)

#gather the fdr p-values from deseq for each timepoint
deseq_res_list = list(res.e, res.t, res.f)
names(deseq_res_list) = c('8hr', '10hr', '14hr')
pull_fdr = function(x){
  data.frame(x) %>% 
    rownames_to_column('gene') %>% 
    dplyr::select(gene, padj)
}
pval_list = map(deseq_res_list, pull_fdr)
map(pval_list, head)
rdat = purrr::reduce(pval_list, full_join, by='gene') %>% 
  column_to_rownames('gene')

#build sig table for limma (0 or 1 for significance)
CUT=0.1
sigdat = apply(rdat < CUT , c(1,2), function(x) as.numeric(x))
sigdat[is.na(sigdat)]<-0

#double-check
apply(sigdat, 2, sum)
sum(res.e$padj < CUT, na.rm=TRUE)
sum(res.t$padj < CUT, na.rm=TRUE)
sum(res.f$padj < CUT, na.rm=TRUE)

#set up venn diagram variables
FILL=c('darkmagenta', 'navy', 'blue')
LABELS=c('8hr', '10hr', '14hr')
df.venn <- data.frame(x = c(-0.866, 0.866, 0),
                      y = c(0.5, 0.5, -1),
                      labels = factor(LABELS, levels=LABELS),
                      fill = as.character(FILL))

#use limma to get venn counts
vdc <- vennCounts(sigdat)
class(vdc) <- 'matrix'
df.vdc <- as.data.frame(vdc)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(-1.5, 1, -0.5, 1, -0.5, 1, 0),
         pcts = paste('(', signif(Counts/sum(Counts) , 3)*100, '%', ')', sep=''),
         labs = paste(Counts, pcts, sep='\n'))


#get solo calls
cdf = df.vdc[,1:3]
solo = df.vdc[apply(cdf,1,sum)==1,]


#build venn
ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .5, size = 0.5, colour = 'black') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom') +
  scale_fill_manual(values = as.character(df.venn$fill)) +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y+0.15, label = df.vdc$Counts, size = 5, color='white', hjust=0.5) +
  annotate("text", x = df.vdc$x, y = df.vdc$y-0.15, label = df.vdc$pcts, size = 3.5, color='white', hjust=0.5) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5), alpha = .5, size = 0.5, colour = 'black')



# format a table for ethanol module genes ---------------------------------

#uncomment to use (loading gene names takes a while)

# ll=load('wgcna/module_membership.Rdata')
# ll
# head(mmdat)
# ll=load('deseq/ethanol_full_LRT_results.Rdata')
# head(res.eth)
# mm = merge(mmdat, data.frame(res.eth), by = 0) %>% 
#   filter(assignment %in% c("darkolivegreen4", "mediumpurple4")) %>% 
#   column_to_rownames('Row.names')
# mm_named = merge_gene_names(mm)
# dim(mm_named)
# mm_named %>% 
#   write_tsv('results/ethanol_module_genes.tsv')

# plot time modules ----------------------------------------------------------------
ll=load('results/timeMods.Rdata')
ll
ll=load('datasets/raw_rld.Rdata')
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

#pu thenm in desired order
rownames(hdat) = hdat$module
ordered_modules = c('MEcoral2',
                   'MEmagenta4',
                   'MEthistle1',
                   'MEmediumpurple4',
                   'MEdarkolivegreen4',
                   'MEhoneydew1')
ohdat = hdat[ordered_modules, ]
ohdat

#plot each of the time module hub genes for ethanol and controls
plotList=list()
for (i in 1:nrow(ohdat)){
  print(paste(i,'...',sep=''))
  row=ohdat[i,]
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


# repeat plotting all other modules ---------------------------------------

modules = colnames(geneModuleMembership)[!colnames(geneModuleMembership) %in% c('gene', timeMods$module)]
geneModuleMembership$gene = rownames(geneModuleMembership)
head(geneModuleMembership)
hubGenes = c()
for (m in modules){
  sub=geneModuleMembership[order(abs(geneModuleMembership[,m])),c(m,'gene')]
  hub=sub[nrow(sub),'gene']
  hubGenes = append(hubGenes, hub)
}
ohdat = data.frame('module'=modules,
                  'hub' = hubGenes,
                  row.names=hubGenes, stringsAsFactors=F) %>% 
  merge_gene_names()
ohdat

#plot each of the time module hub genes for ethanol and controls
plotList=list()
for (i in 1:nrow(ohdat)){
  print(paste(i,'...',sep=''))
  row=ohdat[i,]
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

length(plotList)
plot_grid(plotlist = plotList, nrow=2)





# plot gpc4 ---------------------------------------------------------------

#pick the gene
g = c("ENSDARG00000015472")


#pick the time point
library(DESeq2)
lnames = load('deseq/ethanol_full_LRT_results.Rdata') #All 
resultsNames(dds.eth)
# myPlotCount(dds.eth, g, 'All timepoints')

#run for each individual timepoint
ll = load('deseq/ethanol_8hr_LRT_results.Rdata');print(ll)
ll = load('deseq/ethanol_10hr_LRT_results.Rdata');print(ll)
ll = load('deseq/ethanol_14hr_LRT_results.Rdata');print(ll)

full = my_plotCounts(dds.eth, g, intgroup='treatment', main='All timepoints') + theme(axis.title.y=element_blank())
eight = my_plotCounts(dds.e, g, intgroup='treatment', main='8hr') + theme(axis.title.y=element_blank())
ten = my_plotCounts(dds.t, g, intgroup='treatment', main='10hr') + theme(axis.title.y=element_blank())
fort = my_plotCounts(dds.f, g, intgroup='treatment', main='14hr') + theme(axis.title.y=element_blank())

plot_grid(full, eight, ten, fort, nrow=2, label_x = 'Treatment', label_y='Normalized count', align='v', rel_widths=c(1,1,1,1))


#double-check against original function
plotCounts(dds.eth, g, intgroup='treatment', main='All timepoints')
plotCounts(dds.e, g, intgroup='treatment', main='8hr')
plotCounts(dds.t, g, intgroup='treatment', main='10hr')
plotCounts(dds.f, g, intgroup='treatment', main='14hr')


#check log2 fold changes
res.eth[g,]
res.e[g,]
res.t[g,]
res.f[g,]



# TUKEYS ON FIGURE 7B-------------------------------------------------------------------------

library(readxl)
bdat = read_excel('datasets/my_control_vs_blebbistatin.xlsx')


#get Tukey letters
require(multcomp)
bdat_grp = bdat %>% 
  unite('group', geno, treat) %>% 
  mutate(group = factor(group))
aovgrp <- aov(width ~ group,
              data = bdat_grp)
summary(aovgrp)
tm = glht(aovgrp, linfct = mcp(group = "Tukey"))
tlet = cld(tm)
my_let = tlet$mcletters$monospacedLetters
tukey_df = data.frame(group = names(my_let),
                      tukey = my_let) %>% 
  mutate(tukey = trimws(tukey))

#make plotting df
library(plotrix)
mbdat = bdat %>% 
  group_by(geno, treat) %>% 
  summarize(mn = mean(width),
            se = std.error(width)) %>% 
  ungroup() %>% 
  mutate(treat = factor(treat, levels = c('control', 'blebbistatin')),
         geno = factor(geno, levels = c('wt.wt', 'wt.mut', 'mut.mut')),
         group = paste(geno, treat, sep='_')) %>% 
  left_join(tukey_df, by = 'group') 
  


#build plot
letter_add = 10
mbdat %>% 
  ggplot(aes(x=geno, y=mn, fill=treat)) +
  geom_bar(position='dodge', stat="identity") +
  geom_errorbar(aes(ymin=mn-se, ymax=mn+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  labs(y = 'projection count',
       fill='') +
  scale_fill_manual(values = grey.colors(2)) +
  geom_text(aes(x=geno, y=mn+se+letter_add, label=tukey),
            position=position_dodge(.9))






