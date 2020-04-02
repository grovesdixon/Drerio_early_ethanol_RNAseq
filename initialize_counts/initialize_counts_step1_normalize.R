#DESeq_early_ethanol.R
#Groves Dixon
#8-18-16
#last updated 4-18-19
#This script takes a raw counts table as input
#and outputs the variance stabilized counts,
#uses the WGCNA tutorial methods to identify and remove outliers
#then outputs the variance stabilized counts for job1, job2, 
#and the combined dataset.
#The outputs are used as input for multivariate analyses and WGCNA

#libs
rm(list=ls())
library('DESeq2')




# UPLOAD DATA AND FILTER COUNTS -------------------------------------------------------------


#upload counts table
counts=read.table('./datasets/rnaseq_gene_counts.txt', header=TRUE, row.names=1) #generate this with HTseq (see RNAseq_data_processing_pipeline.txt)

#look at the head of the file
head(counts)
dim(counts)   #32271 features, 55 samples


#how many genes have average count higher than 3?
mns = apply(counts, 1, mean)
table(mns > 3)

#remove non-gene objects from counts data
to.remove=rownames(counts)[grep("__", rownames(counts))]
to.remove
dim(counts)
counts= counts[!rownames(counts) %in% to.remove,]
dim(counts)

#get new totals
tots = apply(counts, 2, sum)
mean(tots)



# BUILD A DATAFRAME ASSOCIATING SAMPLE NAMES WITH TREATMENT CONDIT --------


#get sample.names to parse trait data from
sample.names = colnames(counts)

#set up treatment variable based on sample names
treatment = sample.names
treatment[grep("C", treatment)] <- 'c'
treatment[grep("NoE", treatment)] <- 'c'
treatment[grep("E", treatment)] <- 'e'
table(treatment)  #30 controls and 25 ethanol treated

#set up time variable
time = sample.names
time[grep('6h', time)] <- 6
time[grep('8h', time)] <- 8
time[grep('10h', time)] <- 10
time[grep('14h', time)] <- 14
table(time)   #5 six-hour samples, 20 eight-hour, 20 ten-hour, 10 fourteen-hour

#set up seqjob variable (batches of sequencing)
#*this is a pretty janky way of doing this, but
#should not be a problem when reads are downloaded from SRA with a trait table 
seqjob = sample.names
seqjob[grep('E.10', seqjob)] <- '1' #note the job1 samples have the '.' between E and time
seqjob[grep('E.8', seqjob)] <- '1'  #note the job1 samples have the '.' between E and time
seqjob[grep('NoE', seqjob)] <- '1'  #note job1 samples used 'NoE' to code a control sample
seqjob[grep('C', seqjob)] <- '2'    #job2 controls coded with 'C'
seqjob[grep('E8', seqjob)] <- '2'   #without period
seqjob[grep('E10', seqjob)] <- '2'  #without period
seqjob[grep('E14', seqjob)] <- '2'  #without period
seqjob
table(seqjob) #25 samples in first sequencing job, 30 samples in second sequencing job

#now assemble into a dataframe
coldata <- data.frame(sample.names, treatment, time, seqjob)
#double-check that coldata matches with sample names
coldata
table(coldata$treatment)
table(coldata$time)
write.csv(coldata, file='metadata/coldata.csv', row.names=FALSE)

#write out input data for DESeq scripts
save(counts, coldata, file='./datasets/DEseq_inputs.Rdata')

#------- GET RAW VARIANCE STABILIZED COUNTS ------------

#set up input matrix for DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~1))


#run DESeq
dds = DESeq(ddsHTSeq)

#get DEseq results
res = results(dds)

#get variance stabilized counts and save them
rld = rlog(dds)
rld.df=assay(rld)
colnames(rld.df) = colnames(counts)

#save image of rld dataframe before reducing
save(rld.df, coldata, file='./datasets/raw_rld.Rdata')
save.image(file="datasets/initialize_countsImage1.Rdata")



