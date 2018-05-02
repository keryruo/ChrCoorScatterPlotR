
GenomeFeatureDistribution<-function(origQtlRep,sigfQtlBed,snpBed,geneElementBed,InterGenicBed,plotDir){

        # disable options converting strings to factors
        options(stringsAsFactors = FALSE);

        # FDR<0.1 germline SNP bed file
        commda=paste('awk \'NR>1 { print $0}\'',origQtlRep, '| awk \'($6<0.1) {print $1}\' | awk -F ":|-" \'{print $1"\t"$2"\t"$2+1"\t"$1":"$2}\' | sort -k1,1 -k2,2n > ',sigfQtlBed)
        system(commda) 

        #### Distribution of all germline SNP across genome feature of UTR,exon,intron,intergenic,promoter ####
        # overlap between snp and geneElement
        commdb=paste('bedtools intersect -a', snpBed, '-b', geneElementBed, '-wa -wb | awk \'{ print $8}\' | sort | uniq -c')
        all_num1=read.table(pipe(commdb))
        # overlap between snp and intergenic
        commdc=paste('bedtools intersect -a', snpBed,'-b', InterGenicBed,'-wa -wb | awk \'{ print $8}\' | sort | uniq -c')
        all_num2=read.table(pipe(commdc))
        all_num=rbind(all_num1,all_num2)
        all_num[,1]=all_num[,1]/sum(all_num[,1])

        #### Distribution of germline eQTL at FDR<0.1 level across genome featue of UTR,exon,intron,intergenic,promoter ####
        # overlap between sigfsnp and geneElement
        commdd=paste('bedtools intersect', '-a', sigfQtlBed, '-b', geneElementBed,'-wa -wb | sort | uniq -c | awk \'{ print $9}\' | sort | uniq -c')
        sigf_num1=read.table(pipe(commdd))
        # overlap between sigfsnp and intergenic
        commde=paste('bedtools intersect', '-a', sigfQtlBed, '-b', InterGenicBed, '-wa -wb | sort | uniq -c | awk \'{ print $9}\' | sort | uniq -c')
        sigf_num2=read.table(pipe(commde))
        sigf_num=rbind(sigf_num1,sigf_num2)
        sigf_num[,1]=sigf_num[,1]/sum(sigf_num[,1])

        # construct a df for ggplot
        library(ggplot2)
        df<-rbind(sigf_num,all_num)
        colnames(df)<-c("Proportion","Feature")
        df$Feature<-factor(df$Feature,levels=all_num[order(all_num[,1],decreasing=TRUE),2])

        plotFile=paste0(plotDir,"GenomeFeatureDistribution.tiff")
        tiff(file=plotFile)
        df$Group<-c(rep('eQTLs(FDR<0.1)',nrow(sigf_num)),rep('All SNPs',nrow(all_num)))
        ggplot(df,aes(x=Feature,y=Proportion,fill=Group))+geom_bar(stat='identity',width=0.5,position='dodge')+scale_fill_manual(values=c("#3794bf","#df8640"))  
        ggsave(file=plotFile)
        dev.off()
}



############# Distribution of all-SNPs and eQTLs across peak region #############
InOffPeakDistribution<-function(peakCoorBed,sigfQtlBed,snpBed,plotDir){
        # number of all SNPs locating in peak region
        commda=paste('awk \'{print $0} \'', peakCoorBed ,'| wc -l')
        all_SNP=read.table(pipe(commda))
        commdb=paste('bedtools intersect -a', snpBed,'-b', peakCoorBed,'-wa -wb | awk \'{ print $8}\' | sort | uniq -c | wc -l')
        all_InPeak_SNP=read.table(pipe(commdb))

        # total number of significant eQTLs satisfiying FDR<0.1
        commdc=paste('awk \'{print $0}\'',sigfQtlBed,'| sort | uniq -c| wc -l')
        total_EQTL=read.table(pipe(commdc))
        # intersect peak bed with significant eQTL bed
        commdd=paste('bedtools intersect -a', sigfQtlBed, '-b', peakCoorBed,'-wa -wb | awk \'{{ print $4}}\' | sort | uniq -c | wc -l')
        num_InPeak_EQTL=read.table(pipe(commdd))

        # construct a df for ggplot
        df<-data.frame(Group=c(rep('All SNPs',2),rep('eQTLs(FDR<0.1)',2)),Label=rep(c('in-peak','off-peak'),2),Proportion=c(unlist(all_InPeak_SNP/all_SNP),1-unlist(all_InPeak_SNP/all_SNP),unlist(num_InPeak_EQTL/total_EQTL),1-unlist(num_InPeak_EQTL/total_EQTL)))

        # barplot of germline SNP's distribution across peak
        plotFile=paste0(plotDir,"SigfeQTLPeakDistribution.tiff")
        tiff(file=plotFile)
        ggplot(df, aes(Group,Proportion , fill = Label)) + geom_bar(stat = "identity")  
        ggsave(file=plotFile)
        dev.off()
}



### comparison of gene connectivity between in-peak germEQTL and off-peak germEQTL ###
InOffPeakGeneConnectivity<-function(origQtlRep,peakCoorBed,sigfQtlBed,plotDir){

        # FDR<0.1 SNPs bed file
        commda=paste('awk \'NR>1 { print $0}\'',origQtlRep, '| awk \'($6<0.1) {print $1}\' | awk -F ":|-" \'{print $1"\t"$2"\t"$2+1"\t"$1":"$2}\' | sort -k1,1 -k2,2n > ',sigfQtlBed)
        system(commda) 

        # Distribution of all SNPs
        commdb=paste('bedtools intersect -a', sigfQtlBed, '-b', peakCoorBed, '-wa -wb | awk \'{{ print $4}}\' | sort | uniq')
        InPeak_EQTL=read.table(pipe(commdb))

        commdc=paste('awk \'NR>1 { print $1}\'', origQtlRep, '| awk \'($5<0.1) {print $1}\' |  sort | uniq -c')

        sigfEQTl_AssociatedGeneNum=read.table(pipe(commdc),col.names=c("Associated_gene_num","SNP" ))
        rownames(sigfEQTl_AssociatedGeneNum)<-sigfEQTl_AssociatedGeneNum$SNP
        InPeak_EQTL_AssGeneNum=sigfEQTl_AssociatedGeneNum[unlist(InPeak_EQTL),'Associated_gene_num']
        OffPeak_EQTL_AssGeneNum=sigfEQTl_AssociatedGeneNum[!is.finite(match(rownames(sigfEQTl_AssociatedGeneNum),unlist(InPeak_EQTL))),'Associated_gene_num']

        # boxplot is not appropriate to visualize discrete data,so using barplot instead
        InPeakDf<-data.frame(table(InPeak_EQTL_AssGeneNum))
        colnames(InPeakDf)<-c('genenum','freq')
        InPeakDf$freq<-InPeakDf$freq/sum(InPeakDf$freq)
        OffPeakDf<-data.frame(table(OffPeak_EQTL_AssGeneNum))
        colnames(OffPeakDf)<-c('genenum','freq')
        OffPeakDf$freq<-OffPeakDf$freq/sum(OffPeakDf$freq)
        df<-data.frame(freq=c(InPeakDf$freq,OffPeakDf$freq),genenum=c(InPeakDf$genenum,OffPeakDf$genenum),group=c(rep('in-peak eQTLs',nrow(InPeakDf)),rep('off-peak eQTLs',nrow(OffPeakDf))))
        df$genenum<-factor(df$genenum,levels=0:max(df$genenum))

        plotFile<-paste0(plotDir,"BarPlotGeneConnectivityBetweenInOffPeakGermEQTL.tiff")
        tiff(file=plotFile,width = 1360, height = 720)
        p<-ggplot(df,aes(x=genenum,y=freq,fill=group))+geom_bar(stat='identity',width=0.8,position='dodge')+xlab("Number of genes associated per eQTL")+ylab("Proportion")
        # remove grid and backgroud of the plot region in ggplot2
        p<-p+ theme_bw()+ theme(panel.border = element_blank(),, axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
        # space between tick marks along y-axis
        p<-p+scale_y_continuous(breaks=seq(0,max(df$freq),0.05), limits=c(0,max(df$freq)), expand=c(0,0))
        # legend and other text set
        p+theme(legend.position= "top",axis.title = element_text(family = "Times", color="#666666", face="bold", size=22),axis.text=element_text(size=10),legend.title=element_text(family = "Times", color="#666666", face="bold", size=16),legend.text=element_text(family = "Times", color="#666666", face="bold", size=12)) 
        ggsave(file=plotFile)
        dev.off()
}



#### comparison of sample-sharing number between in-peak germEQTL and off-peak germEQTL ###
InOffPeakEqtlSampleSharing<-function(origQtlRep,snpMatrFile,plotDir){
        library(data.table)
        qtlReport<-fread(origQtlRep,sep="\t",header=TRUE)

        # load all QTLs' matrix
        snpMatr<-read.table(snpMatrFile,header=T,sep="\t",row.names=1)
        # sample-sharing of SNPs
        SNP_sampleShare<-data.frame(luad=apply(snpMatr[,grep("LUAD",colnames(snpMatr))],1,sum),lusc=apply(snpMatr[,grep("LUSC",colnames(snpMatr))],1,sum))
        rownames(SNP_sampleShare)<-colsplit(rownames(snpMatr),"\\-",c('snp','end'))[,'snp']
        # in-peak eQTLs  vs  off-peak eQTLs 
        OffPeak_EQTL<-setdiff(unique(sigfQTLCoor),unlist(InPeak_EQTL))
        InPeak_SampleShar<-SNP_sampleShare[unlist(InPeak_EQTL),]
        OffPeak_SampleShar<-NP_sampleShare[OffPeak_EQTL,]

        luad_sharing<-rbind(data.frame(table(InPeak_SampleShar$luad)),data.frame(table(OffPeak_SampleShar$luad)))
        colnames(luad_sharing)<-c('sample_num','prop')
        luad_sharing$prop<-luad_sharing$prop/sum(luad_sharing$prop)
        luad_sharing$group<-c(rep('in-peak eQTLs',nrow(luad_sharing)/2),rep('off-peak eQTLs',nrow(luad_sharing)/2))

        lusc_sharing<-rbind(data.frame(table(InPeak_SampleShar$lusc)),data.frame(table(OffPeak_SampleShar$lusc)))
        colnames(lusc_sharing)<-c('sample_num','prop')
        lusc_sharing$prop<-lusc_sharing$prop/sum(lusc_sharing$prop)
        lusc_sharing$group<-c(rep('in-peak eQTLs',nrow(lusc_sharing)/2),rep('off-peak eQTLs',nrow(lusc_sharing)/2))

        plotFile<-paste0(plotDir,"LUADSNPSampSharingBetweenInOffPeakGermEQTL.tiff")
        tiff(file=plotFile)
        p<-ggplot(luad_sharing,aes(x=sample_num,y=prop,fill=group))+geom_bar(stat='identity',width=0.8,position='dodge')+xlab("sample number in LUAD")+ylab("Proportion")
        # remove grid and backgroud of the plot region in ggplot2
        p<-p+ theme_bw()+ theme(panel.border = element_blank(),, axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
        # legend and other text set
        p+theme(legend.position= "top",axis.title = element_text(family = "Times", color="#666666", face="bold", size=22),axis.text=element_text(size=10),legend.title=element_text(family = "Times", color="#666666", face="bold", size=16),legend.text=element_text(family = "Times", color="#666666", face="bold", size=12)) 
        ggsave(file=plotFile)
        dev.off()

        plotFile<-paste0(plotDir,"LUSCSNPSampSharingBetweenInOffPeakGermEQTL.tiff")
        tiff(file=plotFile)
        p<-ggplot(lusc_sharing,aes(x=sample_num,y=prop,fill=group))+geom_bar(stat='identity',width=0.8,position='dodge')+xlab("sample number in LUSC")+ylab("Proportion")
        # remove grid and backgroud of the plot region in ggplot2
        p<-p+ theme_bw()+ theme(panel.border = element_blank(),, axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
        # legend and other text set
        p+theme(legend.position= "top",axis.title = element_text(family = "Times", color="#666666", face="bold", size=22),axis.text=element_text(size=10),legend.title=element_text(family = "Times", color="#666666", face="bold", size=16),legend.text=element_text(family = "Times", color="#666666", face="bold", size=12)) 
        ggsave(file=plotFile)
        dev.off()
}

############# genes associated per germline eQTL and germline eQTL associated per gene #############
GeneQTLConnecPieChart<-function(origQtlRep,plotDir){

        ### genes associated per eQTL ###
        # gene connectivity per eQTL #
        commda=paste('awk \'NR>1 { print $0}\'', origQtlRep, '| awk \'($6<0.1) {print $1}\' | sort | uniq -c')
        genenum_perSNP<-read.table(pipe(commda))
        colnames(genenum_perSNP)<-c("num","SNP")
        rownames(genenum_perSNP)<-genenum_perSNP[,2]

        # frequency for genes associated per eQTL
        genenum_prop=c(sum(snpnum_perGene$num==1),sum(genenum_perSNP$num==2), sum(genenum_perSNP$num==3),sum(genenum_perSNP$num>3))/nrow(genenum_perSNP)
        genenum_group=factor(c('1','2','3','>3'),levels=c('1','2','3','>3'))
        genenum_label=paste0(genenum_group,"(",paste0(round(genenum_prop*100),"%",")"))
        .libPaths("/home/tkl/R/x86_64-pc-linux-gnu-library/3.3")
        library(plotrix)
        # piechart plot
        plotFile<-paste0(plotDir,"GeneConnecPerEQTL.tiff")
        tiff(file=plotFile)
        lp1<-pie3D(genenum_prop,labels =genenum_label,border=NA,radius=0.9,col=c("#C8004A","#2B9CD2","#F1E000","#84B61B"),main="Genes connected per germline eQTL")
        dev.off()

        ### eQTL associated per gene ###
        # eQTL connectivity per gene #
        commdb=paste('awk \'NR>1 { print $0}\'', origQtlRep, '| awk \'($6<0.1) {print $2}\' | sort | uniq -c')
        snpnum_perGene<-read.table(pipe(commdb)
        colnames(snpnum_perGene)<-c("num","Gene")
        rownames(snpnum_perGene)<-snpnum_perGene[,2]

        # frequency for  eQTL connectivity per gene
        snpnum_prop=c(sum(snpnum_perGene$num==1),sum(snpnum_perGene$num<=10), sum(snpnum_perGene$num>10 & snpnum_perGene$num<=20),sum(snpnum_perGene$num>20))/nrow(snpnum_perGene)
        snpnum_group=factor(c('1','2-10','11-20','>20'),levels=c('1','2-10','11-20','>20'))
        snpnum_label=paste0(snpnum_group,"(",paste0(round(snpnum_prop*100),"%"),")")
        plotFile<-paste0(plotDir,"eQTLConnecPerGene.tiff")
        tiff(file=plotFile)
        pie3D(snpnum_prop,labels =snpnum_label,border=NA,radius=0.9,labelpos=c(0.7830938,2.3079405,4.9908309,5.9075768),col=c("#D6000E","#0266AC","#E7AA29","#1EA134"),main="germ-eQTLs connected per gene")
        ggsave(file=plotFile)
        dev.off()
}

