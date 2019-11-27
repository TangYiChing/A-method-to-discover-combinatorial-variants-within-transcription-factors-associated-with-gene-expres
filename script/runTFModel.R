suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-s", "--scriptPath", help="path to scripts")
parser$add_argument("-f", "--filePath", help="path to files")
parser$add_argument("-p", "--filePrefix", help="prefix of files")
parser$add_argument("-o", "--outFolder", help="folder name of output")
parser$add_argument("-c", "--chr", type="integer", help="Number of chromosome")
parser$add_argument("-m", "--model", help="type of tfmodel, choices=[tf-binding, tf-regulation, tf-both]")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()
SPATH=args$scriptPath
FPATH=args$filePath
N_CHR=args$chr
FOUT=args$outFolder
PREFIX=args$filePrefix
MODEL=args$model

# use the modified PredictDB_Pipleline_GTEXV6 script
GTEXV6=paste(SPATH, "GTEx_Tissue_Wide_CV_elasticNet_TF_window_zero_index.R", sep="")
source(GTEXV6)
# file path (change this to your own path)
expression_RDS=paste(FPATH, "expression/", PREFIX, ".expression.RDS", sep="")
gene_annot_RDS=paste(FPATH, "annotation/gene_annot/", PREFIX, ".gene_annot.RDS", sep="")
tf_gene_annot_RDS=paste(FPATH, "annotation/tf_gene_annot/", PREFIX, ".tfgene_annot.RDS",sep="")
if (MODEL == 'tf-binding') {
    tf_genotype_dir=paste(FPATH,"genotype/tf-binding_genotype_index/", sep='')
    tf_snp_annot_dir=paste(FPATH, "annotation/tf_snp_annot/nsSNP/",sep="")
} else if (MODEL == 'tf-regulation') {
    tf_genotype_dir=paste(FPATH,"genotype/tf-regulation_genotype_index/", sep='')
    tf_snp_annot_dir=paste(FPATH, "annotation/tf_snp_annot/eQTL/",sep="")
} else if (MODEL == 'tf-both') {
    tf_genotype_dir=paste(FPATH,"genotype/tf-both_genotype_index/", sep='')
    tf_snp_annot_dir=paste(FPATH, "annotation/tf_snp_annot/both/",sep="")
} else { stop("--model should be in [tf-binding, tf-regulation, tf-both]") }
# settings for model
study='tf-model'
n_k_folds=as.numeric('10')
alpha=as.numeric('0.5')
window=as.numeric(1e6)
seed=as.numeric(1345)
logFile <- paste("timelog_chr", N_CHR, ".log", sep="")
# loop through chromosome
for (chr in N_CHR:N_CHR){
  sink(logFile)
  # settings
  chrom=paste(chr,'',sep='')  
  snpset=paste('geuvadis.chr',chr,'snp',sep='')
  out_dir=FOUT

  # get files -- genotype, snp_annot
  genotype_txt=paste(FPATH,'genotype/', PREFIX,'.chr',chr,'.snp.txt', sep='')
  snp_annot_RDS=paste(FPATH, 'annotation/snp_annot/', PREFIX,'.chr',chr,'.snp_annot.RDS', sep='')
   
  # print some progress messages to stderr if "quietly" wasn't requested
  if ( args$verbose ) { 
    cat("Input files:\n")
    cat("    expression file:", expression_RDS, "\n")
    cat("    gene annotation file:", gene_annot_RDS, "\n")
    cat("    genotype file:", genotype_txt, "\n")
    cat("    snp annotation file:", snp_annot_RDS, "\n")
    cat("    tf gene annotation file:", tf_gene_annot_RDS, "\n")
    cat("    tfgenotype directory", tf_genotype_dir, "\n")
    cat("    tfsnp annotation directory", tf_snp_annot_dir, "\n") 
    cat("Model:\n")
    cat("    TF model:", MODEL, "\n")
    cat("    cv=", n_k_folds, "\n")
    cat("    alpha=", alpha, "\n")
    cat("    window (for cis)=", window, "\n")
    cat("    seed number=", seed, "\n")
    cat("Output:\n")
    cat("    output folder:", FOUT, "\n")

  }

  # execute program 
  cat("Runtime\n")
  cat("    chromosome=", N_CHR, "\n")
  start_time <- as.POSIXct(Sys.time())
  TW_CV_model(expression_RDS, genotype_txt, gene_annot_RDS, snp_annot_RDS, tf_gene_annot_RDS,
    n_k_folds, alpha, out_dir, study, chrom, snpset, window, tf_genotype_dir, tf_snp_annot_dir, seed=seed)
  end_time <- as.POSIXct(Sys.time())
  time.taken <- end_time - start_time
  message("Program Start: \t", start_time)
  message("Program End: \t", end_time)
  print(start_time)
  print(end_time)
  print(time.taken)  
}
sink()
