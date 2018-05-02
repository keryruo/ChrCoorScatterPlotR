# a self-defined function to pack common codes
matrixEQTL<-function(snp_file_name,expression_file_name,snp_location_file_name,gene_location_file_name,covariates_file_name,output_file_name_cis,output_file_name_tra){

            library(MatrixEQTL)
            # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
            useModel = modelLINEAR; 
    
            ## Load genotype data
            snp = SlicedData$new();
            snp$fileDelimiter = "\t";      # the TAB character
            snp$fileOmitCharacters = "NA"; # denote missing values;
            snp$fileSkipRows = 1;          # one row of column labels
            snp$fileSkipColumns = 1;       # one column of row labels
            snp$fileSliceSize = 2000;      # read file in slices of 2,000 rows
            snp$LoadFile(snp_file_name);

            ## Load gene expression data
            gene = SlicedData$new();
            gene$fileDelimiter = "\t";      # the TAB character
            gene$fileOmitCharacters = "NA"; # denote missing values;
            gene$fileSkipRows = 1;          # one row of column labels
            gene$fileSkipColumns = 1;       # one column of row labels
            gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
            gene$LoadFile(expression_file_name);
    
            ## load location file of peaks and genes 
            snppos = read.table(snp_location_file_name, header = TRUE, stringsAsFactors = FALSE);
            genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);


            # Only associations significant at this level will be saved
             pvOutputThreshold_cis = 0.005;
             pvOutputThreshold=0.00000001; # trans-eQTL
             output_file_name = output_file_name_tra;
 
             # Error covariance matrix
             # Set to numeric() for identity.
             errorCovariance = numeric();

             # Distance for local gene-peak pairs
             cisDist = 1e6; # 1Mb
             ## Load covariates

             cvrt = SlicedData$new();
             cvrt$fileDelimiter = "\t";      # the TAB character
             cvrt$fileOmitCharacters = "NA"; # denote missing values;
             cvrt$fileSkipRows = 1;          # one row of column labels
             cvrt$fileSkipColumns = 1;       # one column of row labels
             if(length(covariates_file_name)>0) {
             cvrt$LoadFile(covariates_file_name);
             }

             ## Run the analysis ( only perform trans analysis)
             me = Matrix_eQTL_main(snps = snp, gene = gene, cvrt = cvrt,output_file_name = output_file_name,
                      pvOutputThreshold = pvOutputThreshold,useModel = modelLINEAR, 
                      errorCovariance = errorCovariance, verbose = TRUE, 
                      output_file_name.cis = output_file_name_cis,
                      pvOutputThreshold.cis = pvOutputThreshold_cis,
                      snpspos=snppos , genepos = genepos,
                      cisDist = cisDist,pvalue.hist = "qqplot",
                      min.pv.by.genesnp = FALSE,noFDRsaveMemory = FALSE);

             ## Results:
             cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
             cat('Detected local eQTLs:', '\n');
  }
######### parameters setting step ############ 
## passing shell variables to R
# Rscript matrixEQTL.R $snp_file_name $binSignal_file_name $snp_location_file_name $bin_location_file_name $covariates_file_name $output_file_name_cis $output_file_name_tra
args <- commandArgs(trailingOnly = T)
snp_file_name=args[1]
expression_file_name= args[2]
snp_location_file_name=args[3]
gene_location_file_name=args[4]
covariates_file_name=args[5]
output_file_name_cis=args[6]
output_file_name_tra=args[7]

matrixEQTL(snp_file_name,expression_file_name,snp_location_file_name,gene_location_file_name,covariates_file_name,output_file_name_cis,output_file_name_tra)
