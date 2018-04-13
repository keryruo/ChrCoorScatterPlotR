# test data generation
options(stringsAsFactors = FALSE)
set.seed(42)
nchr=22
nsnps=100
test_data2 = data.frame(
       SNP_chr=rep(paste0('chr',1:nchr),each=nsnps),
       SNP_pos=rep(1:nsnps,nchr), 
       gene_chr=rep(paste0('chr',sample(nchr)),nsnps),
       gene_pos=sample(nchr*nsnps)
    )

d<-test_data2
ChrCoorScatterPlot(d)
# pt.col=c("#3794bf", "#df8640")


ChrCoorScatterPlot<-function(d,cex.axis=0.7,pt.col=c('gray10','gray50'),pt.bg=c('gray10','gray50'),pch=21,pt.cex=0.45){
            # Set positions, ticks, and labels for plotting
            d<-d[order(d$SNP_chr, d$SNP_pos),] # sort
            d$SNP_pos2=NA
            d$gene_pos2=NA

            ## Fix for the bug where one chromosome is missing. Adds index column #####
            chrSyb<-as.character(c(1:22,'X','Y','M'))
            # SNP index setting
            # pattern replace,chr10 to 10
            d$SNP_index<-gsub("chr(([0-9]+)|(M)|(X)|(Y))", "\\1", d$SNP_chr)# chr(NO) to NO

            # SNP position setting
            nchr=length(unique(d$SNP_chr))
            SNP_ticks = rep(NA,length(unique(d$SNP_chr))+1)
            SNP_ticks[1] = 0
            for (i in 1:length(unique(d$SNP_index))) {
                 chr=chrSyb[i]
                 d[d$SNP_index==chr, 'SNP_pos2'] <- d[d$SNP_index==chr,'SNP_pos']+SNP_ticks[i]
                 SNP_ticks[i+1] =SNP_ticks[i]+max(d[d$SNP_index==chr,'SNP_pos'])
            }
            xlabel = 'SNP Chromosome'
            # gene index setting 
            # pattern replace,chr10 to 10
            d$gene_index <-gsub("chr(([0-9]+)|(M)|(X)|(Y))", "\\1", d$gene_chr)# chr(NO) to NO

            # gene position setting
            nchr=length(unique(d$gene_chr))
            gene_ticks = rep(NA,length(unique(d$gene_chr))+1)
            gene_ticks[1] = 0
            for (i in 1:length(unique(d$gene_index))) {
                 chr=chrSyb[i]
                 d[d$gene_index==chr, 'gene_pos2'] <- d[d$gene_index==chr,'gene_pos']+gene_ticks[i]
                 gene_ticks[i+1] =gene_ticks[i]+max(d[d$gene_index==chr,'gene_pos'])
            }
            ylabel = 'Gene Chromosome'

            # Initialize plot
            xmax = max(d$SNP_pos2) * 1.03
            xmin = max(d$SNP_pos2) * -0.03
            ymax = max(d$gene_pos2) * 1.03
            ymin = max(d$gene_pos2) * -0.03
            # axis region plot, suppress both xaxis and yaxis
            plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',yaxs='i',xlim=c(xmin,xmax), ylim=c(ymin,ymax),
                 xlab=xlabel,ylab=ylabel,las=1,cex.axis=cex.axis)

            # SNP stagger labels (x-axis)
            # How many chromosome SNP_chr have
            aSNP_chr<-unique(d$SNP_chr)[match(paste0('chr',chrSyb),unique(d$SNP_chr))]
            SNP_labs<-aSNP_chr[!is.na(aSNP_chr)] # existing chr in SNP_chr
            SNP_labs<-c("",SNP_labs)
            SNP_blank = rep('',length(unique(d$SNP_chr))+1)
            # adding tick marks for x-axis
            axis(1,at=SNP_ticks,labels=SNP_blank,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
            axis(1,at=SNP_ticks,labels=SNP_labs,lwd=0,lwd.ticks=0,cex.axis=cex.axis*0.8,line=-0.25)

            # gene stagger labels (y-axis)
            # How many chromosome gene_chr have
            agene_chr<-unique(d$gene_chr)[match(paste0('chr',chrSyb),unique(d$gene_chr))]
            gene_labs<-agene_chr[!is.na(agene_chr)] # existing chr in SNP_chr
            gene_labs<-c("",gene_labs)
            gene_blank = rep('',length(unique(d$gene_chr))+1)
            # adding tick marks for y-axis
            axis(2,at=gene_ticks,labels=gene_blank,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
            axis(2,at=gene_ticks,labels=gene_labs,lwd=0,lwd.ticks=0,cex.axis=cex.axis*0.8)

            # add points per chromosome
            icol=1
            # colour setting
            pt.col = rep(pt.col,length(unique(d$SNP_chr)))[1:length(unique(d$SNP_chr))]
            pt.bg = rep(pt.bg,length(unique(d$SNP_chr)))[1:length(unique(d$SNP_chr))]
            # add all points
            for (i in unique(d$SNP_chr)) {
                 with(d[d$SNP_chr==i, ],points(SNP_pos2, gene_pos2, col=pt.col[icol],bg=pt.bg[icol],cex=pt.cex,pch=pch))
                 icol=icol+1
            }
            # add box for the plot
            box()
}
