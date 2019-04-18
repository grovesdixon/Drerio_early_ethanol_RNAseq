#project_functions.R
library(genefilter)
library(ggplot2)


#modified version of the function provided with the DESeq package
#allows for output of more principal components
mod.plotPCA <- function (object, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n") 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(pca$x[,1:pcs], group = group, 
        intgroup.df, name = colnames(object))
    attr(d, "percentVar") <- percentVar[1:2]
    g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
    y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = 3) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
        100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
        100), "% variance")) + coord_fixed()
   g = g + ggtitle(main)
   g = g + theme_bw()
   print(g)
   return(d)
}



printPCA = function(dat, pc1, pc2){
	ggplot(data = dat, aes_string(x = paste('PC', pc1, sep = ''), y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed()
}




#modified version of the function provided in the DESeq package that can be run from a dataframe
#instead of DESeqTransform object output from
mod.plotPCA.df <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 5, legendTitle='group') 
{
    rv <- rowVars(df)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(df[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    intgroup.df <- as.data.frame(coldat[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        coldat[[intgroup]]
    }
    d <- data.frame(pca$x[,1:pcs], group = group, 
        intgroup.df, name = colnames(df))
    attr(d, "percentVar") <- percentVar[1:2]
    g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
    y = paste('PC', pc2, sep = ''), color = "group")) + 
        geom_point(size = SIZE) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
        100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
        100), "% variance")) + coord_fixed()
   g = g + 
     labs(subtitle=main) +
     guides(color=guide_legend(title=legendTitle))
   # g = g + theme_bw()
   # print(g)
   if (returnData == T){
   return(d)
   } else {
     return(g)
   }
}





write.gomwu.input = function(dat, out.name){
	sign = dat[,'log2FoldChange'] > 0
	sign[sign == TRUE] <- 1
	sign[sign == FALSE] <- -1
	# logp = (-log(dat[,'pvalue'], 10) ) * sign
	stat = dat$stat
	out = data.frame(rownames(dat), stat)
	colnames(out) = c('gene', 'stat')
	head(out)
	write.csv(out, out.name, row.names = F, quote = F)
}




pca.timepoint = function(rld.df, coldata, timepoint){
	eight=rld.df[,colnames(rld.df) %in% coldata$sample.names[coldata$time==timepoint]]
	ec=coldata[coldata$time==timepoint,]
	mod.plotPCA.df(df = eight, coldat = ec, intgroup = 'treatment', main=bquote("Treatment ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)
	mod.plotPCA.df(df = eight, coldat = ec, intgroup = 'seqjob', main=bquote("Batch ("*.(timepoint)*"hr only)"), pc1=1, pc2=2)
	e1=eight[, colnames(eight) %in% ec$sample.names[ec$seqjob==1]]
	e2= eight[, colnames(eight) %in% ec$sample.names[ec$seqjob==2]]
	e1c=ec[ec$seqjob==1,]
	e2c=ec[ec$seqjob==2,]
	e2c$randomized.ethanol=sample(e2c$treatment)
	# mod.plotPCA.df(df = e1, coldat = e1c, intgroup = 'treatment', main=bquote("Treatment (job1 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)#not enough samples because of outlier removal
	mod.plotPCA.df(df = e2, coldat = e2c, intgroup = 'treatment', main=bquote("Treatment (job2 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)
	mod.plotPCA.df(df = e2, coldat = e2c, intgroup = 'randomized.ethanol', main=bquote("Randomized Treatment (job1 "*.(timepoint)*"hr only)"), pc1=1, pc2=2)
}


