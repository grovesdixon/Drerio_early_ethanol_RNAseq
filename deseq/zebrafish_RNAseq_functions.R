#zebrafish_RNAseq_functions.R
library(genefilter)
library(ggplot2)



volcano_plot = function(dat, pcol='pvalue', log2col='log2FoldChange', fdrcol='padj', cut=0.1, XLIM=NULL, YLIM=NULL, MAIN=''){
	plot(-log(dat[,pcol], 10)~dat[,log2col], cex=0.5, xlab=expression(paste("Log"[2], 'Fold Difference')), ylab='-log(p-value)', axes=F, xlim=XLIM, ylim=YLIM, main=MAIN)
	axis(1);axis(2, las=2);box()
	sub=dat[!is.na(dat[,fdrcol]),]
	sub2=sub[sub[,fdrcol]<cut,]
	points(-log(sub2[,pcol], 10)~sub2[,log2col], pch=21,col="black",cex=0.6, bg='red')
}



sig.col='red'
ns.col='black'
addNames = TRUE
topN = 10

ggvolcano_plot = function(deseq.res, sig.col='red', ns.col='black', addNames=F, topN=10, XLIM=F, YLIM=F, MAIN='', submain='', xshift=0.5, yshift=0.5){
	deseq.res=na.omit(data.frame(deseq.res))
	deseq.res$threshold = deseq.res$padj < 0.1
	tolab=deseq.res[1:topN,]
	
	##Construct the plot object
	g = ggplot(data= deseq.res, aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) + scale_colour_manual(values=c(ns.col, sig.col)) +
	  geom_point(alpha=0.4, size=1.75) + xlab("log2 fold difference") + ylab("-log10 p-value") +
	  guides(colour=FALSE) + ggtitle(MAIN, subtitle=submain)
	if (length(XLIM) == 2){
		print("Xlim set")
		g = g + xlim(x=XLIM)
	}
	if (length(YLIM) == 2){
		g = g + ylim(y=YLIM)
	}
	if (addNames){
		print(paste("Gathering gene names for total top genes =", topN))
		sig = merge_gene_names(tolab)
		g=g + geom_text(data=sig, nudge_x=xshift, nudge_y=yshift, aes(x=sig$log2FoldChange, y=-log10(sig$pvalue),
		                     label=sig$external_gene_name), colour="black", size=3) + guides(size=FALSE)
	}
	print(g)
	if (addNames){
		return(sig)
		}
}


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
mod.plotPCA.df <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 5) 
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
   g = g + ggtitle(main)
   g = g + theme_bw()
   print(g)
   if (returnData == T){
   return(d)
   }
}




#Function to write an input file for GO MWU analysis
#dat = a DESeq results object (eg: results(dds, contrast=c('treatment', 'control', 'treated')))
#out.name = the desired name for the output file
write.gomwu.input = function(dat, out.name){
	sign = dat[,'log2FoldChange'] > 0
	sign[sign == TRUE] <- 1
	sign[sign == FALSE] <- -1
	# logp = (-log(dat[,'pvalue'], 10) ) * sign
	stat = abs(dat$stat)*sign
	out = data.frame(rownames(dat), stat)
	colnames(out) = c('gene', 'stat')
	head(out)
	write.csv(out, out.name, row.names = F, quote = F)
}




merge_gene_names = function(dat, col=F, sort.column=F){
	library("biomaRt")
	if(col){
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


get_human_ids = function(geneSet){
	hnames = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description'), filters = c('external_gene_name'), values = geneSet, mart = hsembl)
	hgeneSet=hnames$external_gene_name
	prop=round(length(hgeneSet) / length(geneSet), digits=2)*100
	pct=paste(prop, "%", sep='')
	print(paste("Percentage of gene set with human IDs =", pct))
	return(hgeneSet)
}


get_gene_names = function(geneSet){
	library("biomaRt")
	embl = useMart("ensembl", dataset="drerio_gene_ensembl")
	names = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','description'), filters = c('ensembl_gene_id'), values = geneSet, mart = embl)
	return(names)
}


output_module_for_gomwu = function(dat){
	
}

