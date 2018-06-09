# used for describing significant QTL number under different number of covariates in regression model
# scatter plot

cvrtNoSelect<-function(reaDir,writeDir,N){ 

       library(ggplot2)  
       library(reshape2)
        
       # N is the number of PC
       options(stringsAsFactors = FALSE)
       options(scipen=999)

       stopifnot(class(reaDir)=="character")
       stopifnot(class(writeDir)=="character")

       significantQ<-c()
       for(i in seq(1,N,1)){
           file=paste0(reaDir,"qtlReport_PC",i)
           fileFDR<-paste0(reaDir,"qtlReport_PC",i,"FDR025")
           commd<-paste0('awk \'(NR>1){print $0}\' ',file,' | awk \'($6<0.25){print $0}\' > ',fileFDR)
           system(commd)
           for(fdr in seq(0,0.25,0.01)){
             commd<-paste0('awk \'(NR>1){print $0}\' ',fileFDR,' | awk \'($6<',fdr,'){print $0}\' | wc -l')
             qtlnum<-unlist(read.table(pipe(commd)))
             significantQ<-c(significantQ,qtlnum) 
        }
       }
       sigfPeaks<-data.frame(cvrtnum=factor(rep(paste0(1:N,'PCno'),each=length(seq(0,0.25,0.01))),levels=paste0(1:N,'PCno')),
                                      sigfPeaksNum=significantQ,fdr=rep(seq(0,0.25,0.01),N))
       max_fdr01=max(sigfPeaks[sigfPeaks$fdr==0.1,'sigfPeaksNum'])
       PCno=sigfPeaks[sigfPeaks$sigfPeaksNum==max_fdr01,'cvrtnum']
       PCnum=colsplit(PCno,"PCno",c('num','suffix'))[,'num']
       PCnum=as.numeric(PCnum)
       pdf()
       p0<-ggplot(sigfPeaks,aes(x=fdr,y=sigfPeaksNum,colour=cvrtnum))+geom_line()+
           theme(plot.title = element_text(size=rel(1.5),family="Times",face="bold",hjust=0.5),
           axis.title=element_text(size=9,face="bold",hjust=0.5),axis.text.y=element_text(size=8),
           axis.text.x = element_text(hjust = 1,size=8))
       p0<-p0+ggtitle("Distribution of significant SNP number")+
              xlab("FDR")+ylab("Number of significant SNPs")+ 
              geom_vline(xintercept = 0.1, colour="red", linetype = "dashed")
       p0+geom_text(aes(label=ifelse(sigfPeaksNum==max_fdr01,as.character(cvrtnum),'')),hjust=1,vjust=1,colour='black',size=3)
       ggsave(file=paste0(writeDir,"Significant QTL number vs PC number.pdf"))
       dev.off()
       return(PCnum)
}

# parameters setting
args <- commandArgs(trailingOnly = T)
reaDir=args[1]
writeDir= args[2]
N=as.numeric(args[3])


print(paste("reaDir:",args[1]))
print(paste("writeDir:",args[2]))
print(paste("N:",as.numeric(args[3])))

#function running
cvrtNoSelect(reaDir,writeDir,N)
