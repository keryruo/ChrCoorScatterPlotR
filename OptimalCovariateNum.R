
# used for describing significant QTL number under different number of covariates in regression model
# scatter plot
library(data.table)
library(ggplot2)
cvrtNoSelect<-function(dir,readir){
       significantQ<-c()
       for(i in seq(1,20,1)){
           file=paste0(readir,"qtlReport_PC",i)
           qtlReport<-fread(file,sep="\t",header=TRUE)
           for(fdr in seq(0,0.25,0.01)){
           significantQ<-c(significantQ,sum(qtlReport[,6]<fdr))
        }
       }
       sigfPeaks<-data.frame(cvrtnum=factor(rep(paste0(1:20,'PCno'),each=length(seq(0,0.25,0.01))),levels=paste0(1:20,'PCno')),sigfPeaksNum=significantQ,fdr=rep(seq(0,0.25,0.01),20))
       max_fdr01=max(sigfPeaks[sigfPeaks$fdr==0.1,'sigfPeaksNum'])
       pdf(file=paste0(dir,"Significant QTL number vs PC number.pdf"))
       p0<-ggplot(sigfPeaks,aes(x=fdr,y=sigfPeaksNum,colour=cvrtnum))+geom_line()+theme(plot.title = element_text(size=rel(1.5),family="Times",face="bold",hjust=0.5),axis.title=element_text(size=9,face="bold",hjust=0.5),axis.text.y=element_text(size=8),axis.text.x = element_text(hjust = 1,size=8))
       p0<-p0+ggtitle("Distribution of significant SNP number")+xlab("FDR")+ylab("Number of significant SNPs")+ geom_vline(xintercept = 0.1, colour="red", linetype = "dashed")
       p0+geom_text(aes(label=ifelse(sigfPeaksNum==max_fdr01,as.character(cvrtnum),'')),hjust=1,vjust=1,colour='black',size=3)
       dev.off()
}