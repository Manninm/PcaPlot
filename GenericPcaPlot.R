#This Script will output a PCA for a given gene expression matrix with genes as rows, and samples as columns with an accompanied pheno table, which can consist of <= 4 columns consisting of a directory column (which is used for ballgown), sample-name column, category/disease column, and batch/cohort column. Samples in the Sample Column must match colunmn names in Gene expression matrix. The Script will attempt to clean for strings commonly seen in either HtSeq or Stringtie. Given Batch=FALSE, a pca will be colored by Disease/Condition. If Batch=TRUE, then shapes will be used to identify Batch. The Tissue/Feature options are for easily nameing multiple files via conditions and tissue. 
PcaPlot<-function(ExpMat,transcript=FALSE,batch=FALSE,log=FALSE,Directory,RowNames=1,Skip=1,features='ExperimentalFeatures',tissue='TissueUsedInAnalysis'){
	#ExpMat='TabDelimited Gene Matrix with Genes as rows/Samples as Columns';GroupFile='Phenotype table <= four columns titles Directory, Sample, Disease, Batch';Batch='If true, will use batch column to assign shapes to different batches, must be < 25 batchs';transcript='use if matrix has both geneids and transcript ids assumes txIDs are 2nd column';Directory='use true if your first column consists of directory paths to samples, often used in Ballgown';features='Naming convenience';tissue='Naming Convenience';Log='use if you wish to log2(+1) transform your data.' 
	library(ggfortify)
		exp<-read.table(ExpMat,header=TRUE,row.names=RowNames)
    exp<-exp[,-Skip]
		names(exp)<-gsub("X","",names(exp))
		names(exp)<-gsub('FPKM\\.',"",names(exp))
		names(exp)<-gsub("htseq_nopos.inter.str.txt","",names(exp))
		names(exp)<-gsub("\\_htseqout_noposuniontxt","",names(exp))
		names(exp)<-gsub('FPKM.Sample_','',names(exp))
		names(exp)<-gsub('Sample_','',names(exp))
		names(exp)<-gsub('*_merged',"",names(exp))
		names(exp)<-gsub("\\.","-",names(exp))
		names(exp)<-gsub("X","",names(exp))
		 exp[is.na(exp)] <- 0 
		 exp.pca<-prcomp(t(exp))
		 na<-names(exp)
		 pdf(paste("PCA1v2&2v3&3v4_",tissue,features,".pdf",sep=""), width=15, height=15)
		 plot(exp.pca$x[,1], exp.pca$x[,2], main=paste("PCA1vs2SRC_", features, tissue,sep=""), xlab = "PCA 1", ylab = "PCA 2")
		 text(exp.pca$x[,1], exp.pca$x[,2], labels=na, pos= 3) #labels points
		 plot(exp.pca$x[,2], exp.pca$x[,3], main=paste("PCA2vs3_SRC", features, tissue,sep=""), xlab = "PCA 2", ylab = "PCA 3")
		 text(exp.pca$x[,2], exp.pca$x[,3], labels=na, pos= 3) #labels points
		 plot(exp.pca$x[,3], exp.pca$x[,4], main=paste("PCA3vs4_SRC", features, tissue,sep=""), xlab = "PCA 3", ylab = "PCA 4")
		 text(exp.pca$x[,3], exp.pca$x[,4], labels=na, pos= 3) #labels points
		 dev.off()
}