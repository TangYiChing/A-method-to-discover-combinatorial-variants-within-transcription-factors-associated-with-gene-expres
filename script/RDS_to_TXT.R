#!/usr/bin/R

## This is designed for converting expression file 
# .RDS: sample_id by gene_id
# .txt: gene_id by sample_id



#Convert RDS to txt format or the other way round
#
#    file: expression file
#    method: 1: yes, RDS to txt; 0: no, txt to RDS


# get args
argv <- commandArgs(trailingOnly = TRUE)
expressionfile <- argv[1]
outFile <- argv[2]
method <- as.integer(argv[3])


# main
if (method == 1){
    #message("method = 1, RDS to txt") 
    # Load data
    expression <- readRDS(expressionfile)
    # Transpose expreesion
    #expression <- t(expression)
    # Save to txt
    write.table(expression, file=outFile, sep="\t", quote=F,  row.names=FALSE)

}else if (method == 0){
    #message("method = 0, txt to RDS")
    # Load data
    expression <- read.table(expressionfile, sep="\t", stringsAsFactors = FALSE, header = TRUE)
    # Transpose expression.
    #expression <- t(expression)
    # Save to RDS
    saveRDS(expression, outFile)

} else { 
      stopifnot("method should be 0 or 1, given=", method)
    
}

#RDSout <- argv[2]
#expressionfile_out <- argv[3]

