######' wgcna gene module and trait association plot ######'
######' including go enrichment barplot, boxplot of FPKM difference for goups of module-associated trait ######'
######' 2018/4/26 ######'


## ##' ggplot theme setting
## ##' @title 
## ##' @param font.size font size

## theme of barplot displaying GO enrichment 
theme_barplot <- function(font.size=14) {
  theme_bw() +
    theme(axis.text.x = element_text(colour = "#666666",family="Times",face="bold",
                                     size = font.size+1, vjust =1 ),
          axis.text.y = element_text(colour ="#666666",
                                     size = font.size, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color ="#666666",face="bold",
                                    size = font.size+2),
          axis.title.y = element_text(angle=90,family="Times",face="bold",size=font.size+1)
    )
}

## theme of boxplot displaying FPKM difference in groups of module-associated trait                               
theme_boxplot <- function(font.size=12) {
  theme_bw() +        
          theme(plot.title = element_text(family = "Times", color="#666666", face="bold",size=font.size+6,hjust=0.5),
           axis.title = element_text(family = "Times", color="#666666", face="bold", size=font.size+3),
           axis.text.x=element_text(family = "Times", color="#666666", face="bold", size=font.size) ,
           legend.title=element_text(family = "Times", color="#666666", face="bold", size=font.size+3),
           legend.text=element_text(family = "Times", color="#666666", face="bold", size=font.size))
}


######### load packages #######
library(pheatmap)
library(dendsort)
library(ggplot2)
library(gridExtra)

##################  input description ################ 
#### module_trait: a list with gene module as names and associated traits as elements
#### moduleDavidDir: module DAVID annotation file
#### datExpr: gene expression matrix with gene in columns and sample in rows

moduleTraitAssociationPlot<-function(module_trait,moduleDavidDir,datExpr){
  for(i in 1:length(module_trait)) {
    plot_list=list()
    # module name
    module<-names(module_trait)[i]
    module_genes<-rownames(geneModule)[is.finite(match(geneModule$col,module))]
    submatr<-datExpr[,module_genes]
    boxdata<-melt(submatr)
    colnames(boxdata)<-c("sample","gene","FPKM")
    boxdata$logFPKM<-log(boxdata$FPKM+1,2)
    
    # module pheatmap of FPKM
    pdf()
    h<-pheatmap(t(submatr), color =colorRampPalette("#009999","#FF6666")(100), 
              annotation_col=cliPhe33,annotation_colors =ann_colors,scale='row',
              cluster_rows= as.hclust(dendsort(hclust(dist(t(submatr))))),cluster_cols = F,
              show_rownames = F,main=paste0("MM",module))
    dev.off()
    
    ## module DAVID annotation barplot
    moduleFile<-paste0(moduleDavidDir,module,"_DAVID.txt")
    enrich<-read.table(moduleFile,header=T,sep="\t")
    enrich<-enrich[enrich$FDR<0.2,c('Category','Term','PValue',"Fold.Enrichment","Bonferroni","FDR" )]
    enrich$logP<--log(enrich$PValue,10)
    enrich<-enrich[order(enrich$logP),]
    enrich$Term<-factor(enrich$Term,levels=enrich$Term)
    
    p <- ggplot(enrich, aes(x = Term, y = logP))+ylab("-logP")
    p <- p + geom_bar(stat = "identity") + coord_flip() + theme_barplot(12)
    p <- p + aes(fill=logP) + scale_fill_continuous(low="#009999", high="#FF6666")
    
    plot_list[[1]]=h[[4]]
    plot_list[[2]]=p
    
    for(j in 1:length(module_trait[[i]])){
        trait<-module_trait[[i]][j]
        boxdata[,trait]<-cliPhe33[boxdata$sample,trait]
        boxdata[,trait]<-as.factor(boxdata[,trait])
        
        ## boxplot of FPKM in different groups of module-associated trait
        p<-ggplot(boxdata,aes(x=boxdata[,trait],y=logFPKM,fill=boxdata[,trait]))+geom_boxplot()
        p<-p+scale_color_manual(values =c("#009999","#FF6666"))
        p<-p+ggtitle(paste0("MM",module,"_",trait))+xlab(trait)+guides(fill=guide_legend(title=trait))
       
        # remove grid and backgroud of the plot region in ggplot2
        p<-p+ theme_bw()+ theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), 
           panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
        p<-p+theme_boxplot(12)
        plot_list[[j+2]]=p
    }
    pdf()
    g<-do.call(grid.arrange,plot_list)
    ggsave(g,file=paste0(heatmapDir,names(module_trait)[i],".pdf"),width=560, height=560,units = "mm")
    dev.off()
  }
}
      
