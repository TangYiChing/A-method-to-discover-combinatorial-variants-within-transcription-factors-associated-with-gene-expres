#!/usr/bin/Rscript
###
###
#library.path <- .libPaths()
#library("timeseries", lib.loc = library.path)
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("car")) # outlierTest, infindexPlot
suppressPackageStartupMessages(library("moments")) # kurtosis
suppressPackageStartupMessages(library("ggpubr")) # multiple plots
suppressPackageStartupMessages(library("gridExtra")) # aarangeGrob
suppressPackageStartupMessages(library("ggfortify"))



read_as_df <- function(y_table_fname){
  #.libPaths("/home/ytang4/anaconda3/lib/R/library")

  df <- read.table(y_table_fname, sep="\t", header=TRUE, stringsAsFactors=F)
  df
}

plot_r2 <- function(df, significant_outliers){
  #.libPaths("/home/ytang4/anaconda3/lib/R/library")


  #initialize
  fit <- lm(formula = TF_R2 ~ baseline_R2, data=df)
  sub_df = df[which(df$TF_R2 > df$baseline_R2),] # only interested in datapoints where TF>PrediXcan
  n_size = nrow(sub_df)

  #plot
  r2 = summary(fit)$r.squared
  r2 = round(r2,3)
  scatter.outliers <- ggplot(df, aes(x=baseline_R2, y=TF_R2)) +
                     geom_point(shape=1) +
                     geom_smooth(method=lm, se=FALSE) +
                     theme_bw() +

                     ggtitle(paste(chrom, "\n", "R-squared=", r2, "(n=", n_size, ")")) +
                     geom_point(data=sub_df[significant_outliers,], aes(x=baseline_R2, y=TF_R2), size=0.5, colour="red") +
                     theme(plot.title = element_text(hjust = 0.5, size=12))
                     #geom_text(data=sub_df[significant_outliers,], aes(x=baseline_R2, y=TF_R2), label=sub_df[significant_outliers,]$genename, color = "red", vjust = -0.5, size = 3) + theme(plot.title = element_text(hjust = 0.5, size=12))
  figure = scatter.outliers
  outf = paste(chrom, "outlierPlot.png", sep=".")
  ggsave(filename=outf, plot=figure, width = 4, height = 4, dpi = 150, units = "in", device='png')
  
  message("Outlier Plot= ", outf)
}

get_outlier <- function(df){
  #.libPaths("/home/ytang4/anaconda3/lib/R/library")


  # fit a linear model
  fit <- lm(formula = TF_R2 ~ baseline_R2, data=df)
  # calculate studentized residuals
  sresid <- rstudent(fit)
  df.sresid <- as.data.frame(sresid)
  # distribution of residuals (with normality check using Kurtosis)
  #kvalue <- kurtosis(df.sresid$sresid)
  # detect outliers (significant outliers are defined as |residuals| > 2 and Bonferroni P value < 0.01)
  p_values <- outlierTest(fit, cutoff=0.01, n.max=nrow(df), labels=names(rstudent(fit))) #cutoff=Inf and n.max=Inf to show all observations. [Note] original setting: cutoff=0.01, n.max=nrow(df)
  significant_outliers <- names(p_values$rstudent) # index number, e.g. 44, 41
  rstudent <- p_values$rstudent
  unadjusted.p_value <- p_values$p
  Bonferonni.p_value <- p_values$bonf.p
  outliers.df <- data.frame(rstudent, unadjusted.p_value, Bonferonni.p_value)
  # combine dataframe for output
  sub_df = df[which(df$TF_R2 > df$baseline_R2),] # only interested in datapoints where TF>PrediXcan
  new <- cbind(sub_df[significant_outliers,], outliers.df[significant_outliers,])
  new.omna <- na.omit(new)
  n_outliers <- nrow(new.omna)
  outlier_f = paste(chrom, n_outliers,"outliers.txt", sep=".")
  write.table(new.omna, file=outlier_f, quote=FALSE, sep="\t", row.names=FALSE)
  message("Outlier Table=", outlier_f)
  significant_outliers
}

main <- function(ytable, out_prefix){

  #.libPaths("/home/ytang4/anaconda3/lib/R/library")

  df = read_as_df(ytable)
  sig_outliers  = get_outlier(df)
  plot_r2(df, sig_outliers)
}


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("--ytable", type="character", 
    help="File name of ytable. Headers=[gene, baseline_R2, TF_R2, genename]")
parser$add_argument("--out_prefix", type="character", default="result", 
    help="Output prefix")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

y_table = args$ytable
chrom = args$out_prefix
main(y_table, chrom)
