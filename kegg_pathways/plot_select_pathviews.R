# #plot_select_pathviews.R
# 
# 
# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()
library(dplyr)
library(DESeq2)
library(tidyverse)

ll=load('kegg_pathways/pathview_input.Rdata') #build this with kegg_tests.R

#cherry-picked set
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
setwd('./kegg_pathways/selected_pathview_figures/')
tmp = sapply(keggresids,
             function(pid) 
               pathview(gene.data=foldchanges,
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
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="dre", node.sum= NODE.SUM, low= LOW, mid=MID, high=HIGH, limit=list(gene=scale.limit,cpd=scale.limit), kegg.native=TRUE))
setwd('../../')


# assemble keggs and genes for checking -----------------------------------

#assemble a single dataframe
kegg_vector =  unlist(kegg.sets.dr[keggresids])
gene_df = keggs %>% 
  dplyr::rename(entrez = gene) %>% 
  filter(kegg %in% keggresids) %>% 
  dplyr::inner_join(data.frame(res), by = 'entrez') %>% 
  as_tibble()


# look for a particular kegg --------------------------------------------------

#kegg subset
my_kegg = 'dre04310'
kdat = gene_df %>% 
  dplyr::filter(kegg==my_kegg)
kdat


# write out dataframes for each selected kegg -----------------------------

for (i in 1:length(keggresids)){
  my_kegg = keggresids[i]
  kegg_name = knames[i]
  print(paste(kegg_name, '...', sep=''))
  out_name = paste(my_kegg, kegg_name, 'gene_results.csv', sep='_')
  out_path = paste('kegg_pathways/selected_pathview_figures/', out_name, sep='')
  kdat = gene_df %>% 
    dplyr::filter(kegg==my_kegg)
  kdat %>% 
    write_csv(out_path)
}


