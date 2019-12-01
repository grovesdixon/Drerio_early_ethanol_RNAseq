#get_Drerio_gene_names.R

#upload a table with ensemble ids you want names for
dataToName = read.table('~/gitreps/vangl2_gene_expression/deseq/ethanol_results.tsv')

#set the column name to ensembl_gene_id
dataToName$ensembl_gene_id = rownames(dat)


#This is the function that runs biomaRt
merge_gene_names = function(dat, col=F, sort.column=F){
  library("biomaRt")
  if(col!=FALSE){
    geneSet=dat[,col]
  }
  else{
    geneSet=rownames(dat)
  }
  embl = useMart("ensembl", dataset="drerio_gene_ensembl")
  names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description', 'chromosome_name', 'start_position', 'end_position'), filters = c('ensembl_gene_id'), values = geneSet, mart = embl)
  dat2=dat
  dat2$ensembl_gene_id<-geneSet
  dat3=merge(dat2, names, by = 'ensembl_gene_id', all.x=T)
  if(sort.column!=F){
    dat3=dat3[order(dat3[,sort.column]),]
  }
  return(dat3)
}



#Run the funciton on the dataset
ndat = merge_gene_names(dat = dataToName)