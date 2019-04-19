#venn_diagram_functions.R
library(VennDiagram)

get.sig=function(res, TIME){
  sig=data.frame(res[!is.na(res$padj) & res$padj<=CUT,])
  sig$time=TIME
  x=nrow(sig)
  print(paste(x, 'significant genes'))
  sig=sig[order(sig$pvalue),]
  return(sig)
}

pairwiseVenn=function(dat1, dat2, lab1, lab2){
	grid.newpage()
	draw.pairwise.venn(area1=nrow(dat1), area2=nrow(dat2), cross.area=length(intersect(rownames(dat1), rownames(dat2))), category = c(lab1, lab2), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = T, euler.d=T)
}


threeway_venn = function(s1, s2, s3, CATEGORY){
  i12 = intersect(s1,s2)
  n12=length(intersect(s1,s2))#12
  i13=intersect(s1, s3)
  n13=length(intersect(s1, s3))#13
  i14=intersect(s1, s4)
  n14=length(intersect(s1, s4))#14
  i23=intersect(s2, s3)
  n23=length(intersect(s2, s3))#23
  i24=intersect(s2, s4)
  n24=length(intersect(s2, s4))#23
  i123=intersect(  intersect(s1, s2), s3)
  n123=length(intersect(  intersect(s1, s2), s3) )#123
  grid.newpage()
  draw.triple.venn( area1 = nrow(sig.e),
                  area2 = nrow(sig.t),
                  area3 = nrow(sig.f),
                  n12=n12,
                  n13=n13,
                  n14=n14,
                  n23=n23,
                  n24=n24,
                  n123=n123,
                  category=CATEGORY,
                  fill=c("skyblue", "pink1", "mediumorchid")
  )
  res=list(i12, i13, i23, i123)
  names(res) = c('i12', 'i13', 'i23', 'i123')
  return(res)
}


fourway_venn = function(s1, s2, s3, s4, CATEGORY){
	i12 = intersect(s1,s2)
	n12=length(intersect(s1,s2))#12
	i13=intersect(s1, s3)
	n13=length(intersect(s1, s3))#13
	i14=intersect(s1, s4)
	n14=length(intersect(s1, s4))#14
	i23=intersect(s2, s3)
	n23=length(intersect(s2, s3))#23
	i24=intersect(s2, s4)
	n24=length(intersect(s2, s4))#23
	i34=intersect(s3, s4)
	n34=length(intersect(s3, s4))#24
	i123=intersect(  intersect(s1, s2), s3)
	n123=length(intersect(  intersect(s1, s2), s3) )#123
	i124=intersect(  intersect(s1, s2), s4)
	n124=length(intersect(  intersect(s1, s2), s4))#124
	i134=intersect(  intersect(s1, s3), s4)
	n134=length(intersect(  intersect(s1, s3), s4))#134
	i234=intersect(  intersect(s2, s3), s4)
	n234=length(intersect(  intersect(s2, s3), s4))#234
	i1234=intersect( intersect(  intersect(s1, s2), s3), s4)
	n1234=length(intersect( intersect(  intersect(s1, s2), s3), s4))#1234
	grid.newpage()
	draw.quad.venn( area1 = nrow(sig.all),
		area2 = nrow(sig.e),
		area3 = nrow(sig.t),
		area4 = nrow(sig.f),
		n12=n12,
		n13=n13,
		n14=n14,
		n23=n23,
		n24=n24,
		n34=n34,
		n123=n123,
		n124=n124,
		n134=n134,
		n234=n234,
		n1234=n1234,
		category=CATEGORY,
		fill=c("skyblue", "pink1", "mediumorchid", "orange")
	)
	res=list(i12, i13, i14, i23, i34, i123, i124, i134, i234, i1234)
	names(res) = c('i12', 'i13', 'i14', 'i23', 'i34', 'i123', 'i124', 'i134', 'i234', 'i1234')
	return(res)
}