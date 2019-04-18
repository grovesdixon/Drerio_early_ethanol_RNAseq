#figure_plotting.R
#plot the figures for publication

library(tidyverse)
library(cowplot)

# SUMMARY PLOTS -----------------------------------------------------------
#heatmap and pcas

lnames=load("datasets/large_ignored/initialize_countsImage1.Rdata")
source("./initialize_counts/multivariate_functions.R")

#### MAIN HEATMAP
x=colnames(rld.df)
y=sub("NoE", "C", x)
z=sub("h.", "", y)
l = substr(z, 1,nchar(z)-1)
l2=sub(".", "", l, fixed=T)
pheatmap(cor(rld.df, method = 'spearman'), labels_row=l2, labels_col=labs, treeheight_row=0, treeheight_col=0, number_color='blue')


#### PCAS
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

