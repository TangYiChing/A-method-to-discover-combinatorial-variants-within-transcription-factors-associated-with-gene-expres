suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-s", "--scriptPath", help="path to scripts")
parser$add_argument("-f", "--filePath", help="path to files")
parser$add_argument("-o", "--outFolder", help="folder name of output")
parser$add_argument("-c", "--chr", type="integer", help="Number of chromosome")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
SPATH=args$scriptPath
FPATH=args$filePath
N_CHR=args$chr
FOUT=args$outFolder

# use the modified PredictDB_Pipleline_GTEXV6 script
GTEXV6=paste(SPATH, "GTEx_Tissue_Wide_CV_elasticNet.R", sep="")
source(GTEXV6)
# get files
#FPATH='/data2/GTEx_v6/code/script/TFs_model_w_GTEx_Tissue_Wide_CV_elasticNet/geuvadis_data/'
expression_RDS=paste(FPATH, "expression/geuvadis.expression.RDS", sep="")
gene_annot_RDS=paste(FPATH, "annotation/gene_annot/geuvadis.gene_annot.RDS", sep="")
# settings for model
study='baseline-model'
n_k_folds=as.numeric('10')
alpha=as.numeric('0.5')
window=as.numeric(1e6)
seed=as.numeric(1345)
logFile <- paste("timelog_chr", N_CHR, ".log", sep="")
# loop through chromosome
for (chr in N_CHR:N_CHR){
  sink(logFile)
  print(paste("The is chromosome", chr))
  # settings
  chrom=paste(chr,'',sep='')  
  snpset=paste('geuvadis.chr',chr,'snp',sep='')
  out_dir=FOUT

  # get files -- genotype, snp_annot
  genotype_txt=paste(FPATH,'genotype/geuvadis.chr',chr,'.snp.txt', sep='')
  snp_annot_RDS=paste(FPATH, 'annotation/snp_annot/geuvadis.chr',chr,'.snp_annot.RDS', sep='')
   
  # execute program 
  start_time <- as.POSIXct(Sys.time())
  TW_CV_model(expression_RDS, genotype_txt, gene_annot_RDS, snp_annot_RDS,
    n_k_folds, alpha, out_dir, study, chrom, snpset, window, seed=seed)
  end_time <- as.POSIXct(Sys.time())
  time.taken <- end_time - start_time
  message("Program Start: \t", start_time)
  message("Program End: \t", end_time)
  print(start_time)
  print(end_time)
  print(time.taken)  
}
sink()
