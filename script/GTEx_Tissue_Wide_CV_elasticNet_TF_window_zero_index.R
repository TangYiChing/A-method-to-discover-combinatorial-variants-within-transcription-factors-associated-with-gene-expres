suppressMessages(library(glmnet))
suppressMessages(library(methods))
"%&%" <- function(a,b) paste(a, b, sep = "")

TW_CV_model <- function(expression_RDS, genotype_txt, gene_annot_RDS, snp_annot_RDS, tf_gene_annot_RDS, n_k_folds, alpha, out_dir, tis, chrom, snpset, window, tf_genotype_dir, tf_snp_annot_dir, seed = NA) {
  # Load data
  expression <- readRDS(expression_RDS)
  class(expression) <- 'numeric'
  genotype <- read.table(genotype_txt, header = TRUE, row.names = 'Id', stringsAsFactors = FALSE)
  genotype <- t(genotype)  #Transpose genotype for glmnet
  gene_annot <- readRDS(gene_annot_RDS)
  gene_annot <- subset(gene_annot, gene_annot$chr == chrom)
  snp_annot <- readRDS(snp_annot_RDS)
  tf_annot <- readRDS(tf_gene_annot_RDS)
  tf_annot <- subset(tf_annot, tf_annot$gene_chr == chrom)
  rownames(gene_annot) <- gene_annot$gene_id
  # Subset expression data to only include genes with gene_info
  sharedGenes <- intersect(colnames(expression), rownames(gene_annot))
  if(length(sharedGenes) > 1){
      expression <- expression[,sharedGenes]
  }else {
      expression <- subset(expression, select=c(sharedGenes))
  }
  #expression <- expression[, intersect(colnames(expression), rownames(gene_annot))]
  exp_samples <- rownames(expression)
  exp_genes <- colnames(expression)
  n_samples <- length(exp_samples)
  n_genes <- length(exp_genes)
  seed <- ifelse(is.na(seed), sample(1:2016, 1), seed)
  # Output files
  log_df <- data.frame(chrom, n_genes, seed, alpha)
  colnames(log_df) <- c('chr', 'n_genes', 'seed_for_cv', 'alpha')
  write.table(log_df, file = out_dir %&% tis %&% '_chr' %&% chrom %&% '_elasticNet_model_log.txt', quote = FALSE, row.names = FALSE, sep = "\t")
  set.seed(seed)
  groupid <- sample(1:n_k_folds, length(exp_samples), replace = TRUE)

  resultsarray <- array(0, c(length(exp_genes), 13))  # add w_tf, pct_tf, ttl_snps, used_tf
  if(length(exp_genes)==1){
      dimnames(resultsarray)[[1]] <- list(exp_genes)}
  else{
      dimnames(resultsarray)[[1]] <- exp_genes
  }
  resultscol <- c("gene", "alpha", "cvm", "lambda.iteration", "lambda.min", "n.snps", "R2", "pval", "genename", "w_TF", "pct_TF", "total.snps", "used_TF")
  dimnames(resultsarray)[[2]] <- resultscol
  workingbest <- out_dir %&% "working_TW_" %&% tis %&% "_exp_" %&% n_k_folds %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% ".txt"
  write(resultscol, file = workingbest, ncolumns = 13, sep = "\t")

  weightcol <- c("gene","varID", "rsid","ref","alt","beta","alpha")
  workingweight <- out_dir %&% "TW_" %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% ".txt"
  write(weightcol, file = workingweight, ncol = 7, sep = "\t")

  covariance_out <- out_dir %&% tis %&% '_chr' %&% chrom %&% '_snpset_' %&% snpset %&% '_alpha_' %&% alpha %&% "_covariances.txt"
  
  # Training gene by gene
  for (i in 1:length(exp_genes)) {
    gene <- exp_genes[i]
    #message(i,"/",length(exp_genes), "gene=", gene)
    # Reduce genotype data to only include SNPs within specified window of gene.
    geneinfo <- gene_annot[gene,]
    start <- geneinfo$start - window
    end <- geneinfo$end + window
    # Pull cis-SNP info
    cissnps <- subset(snp_annot, snp_annot$pos >= start & snp_annot$pos <= end)
    # Pull cis-SNP genotypes
    cisgenos <- genotype[,intersect(colnames(genotype), cissnps$varID), drop = FALSE]
    num_ori_cissnps <- ncol(cisgenos)
    # Search for associated transcription factors
    tf_info <- subset(tf_annot, gene_id==gene) # for transsnps
    n_tf <- nrow(tf_info)
    message(i,"/",length(exp_genes), ", gene=", gene, ", #TF=", n_tf)
    if( ncol(cisgenos)>0 & n_tf>0 ){
        w_tf = "yes"
        num_transsnps <- 0
        for( j in 1:n_tf){  # loop through each tfs for the same gene
            fname_tf_geno=Sys.glob(file.path(tf_genotype_dir, paste0(tf_info[j,]$TF_name, '.', 'chr', tf_info[j,]$TF_chr, '.', '*'))) # delimeter of filename = '.'
            fname_tf_snp_annot=Sys.glob(file.path(tf_snp_annot_dir, paste0('*', 'chr', tf_info[j,]$TF_chr, '.', '*')))
            message(j, "/", n_tf, "TFs", " name=", tf_info[j,]$TF_name) #", fname_tf_genotype file=", fname_tf_geno, ", fname_tf_snp_annot file=", fname_tf_snp_annot)
            #if(isTRUE(!file.exists(fname_tf_geno))){
            #    next
            #    }
            if(isTRUE(file.exists(fname_tf_geno))){
                used_tf = paste(tf_info$genetfpair, sep=" ", collapse=",")
                # found genotype and snp_annot files and load data
                tf_geno <- read.table(fname_tf_geno, sep="\t", header = TRUE, row.names = 'Id', stringsAsFactors = FALSE)
                transgenos <- t(tf_geno)  #Transpose genotype for glmnet
                tf_snp_annot <- read.table(fname_tf_snp_annot, sep="\t", header = TRUE, stringsAsFactors = FALSE)
                rownames(tf_snp_annot) <- tf_snp_annot$varID
                
                transsnps <- tf_snp_annot[colnames(transgenos), ]#subset(tf_snp_annot, tf_snp_annot$pos >= tf_start & tf_snp_annot$pos <= tf_end)
                snp_annot <- rbind(snp_annot, transsnps) # cat transsnps to snp_annot for betas annotation
                snp_annot <- snp_annot[!duplicated(rownames(snp_annot)),]
                cisgenos <- cbind(cisgenos, transgenos)
                cisgenos <- cisgenos[, !duplicated(colnames(cisgenos))] # remove duplicate snps
                num_transsnps <- num_transsnps + ncol(transgenos) 
            }else{
                message("    tf= ", tf_info[j,]$TF_name,"has no tf snps") # has TFs but NO snps
            }
        } # end-of-adding TFs
    }
    else{
        w_tf = "no"
        num_transsnps <- 0
        used_tf = "None"
    }
    
    # Reduce cisgenos to only include SNPs with at least 1 minor allele in dataset
    #message(dim(cisgenos))
    #cm <- colMeans(cisgenos, na.rm = TRUE)
    #message(dim(colMeans(cisgenos)), head(colMeans(cisgenos)))
    #minorsnps <- subset(colMeans(cisgenos), cm > 0 & cm < 2)
    #minorsnps <- names(minorsnps)
    #cisgenos <- cisgenos[,minorsnps, drop = FALSE]
    maf<- as.numeric(0.01)
    effect_allele_freqs <- colMeans(cisgenos[ , drop = FALSE]) / 2
    cisgenos <- cisgenos[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf))]
    ttl_snps <- ncol(cisgenos) #total candidate snps


    model.start <- as.POSIXct(Sys.time())
    if ( is.null(ncol(cisgenos)) || ncol(cisgenos) < 2 ) {
      #message(ncol(cisgenos))
      # Need 2 or more cis-snps for glmnet.
      bestbetas <- data.frame()
      pct_tf <- 0
    } else {
      # Pull expression data for gene
      exppheno <- expression[,gene]
      # Scale for fastLmPure to work properly
      exppheno <- scale(exppheno, center = TRUE, scale = TRUE) 
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(expression)      
      # Run Cross-Validation
      # parallel = TRUE is slower on tarbell
      #message("        #current gene: ", gene, '\t', i, "/", length(exp_genes))
      #message("        #has transsnps (TF addition)=\t", w_tf)
      #message("        #cross-validation...")
      #message("        #size of expression\t", '\t#samples=' ,nrow(exppheno), '\t#genes=', ncol(exppheno))
      #message("        #size of cisgenos\t",  '\t#samples=' ,nrow(cisgenos), '\t#snps=', ncol(cisgenos))
      #message("        #size of transgenos\t", num_transsnps, "\t(",num_transsnps/ncol(cisgenos)*100, "%)")
      #message("        #used gene-tf\t", used_tf)
      pct_tf <- num_transsnps/ncol(cisgenos) # percentage of tf in candidate genotypes for modeling
      bestbetas <- tryCatch(
        { fit <- cv.glmnet(as.matrix(cisgenos),
                          as.vector(exppheno),
                          nfolds = n_k_folds,
                          alpha = alpha,
                          keep = TRUE,
                          foldid = groupid,
                          parallel = FALSE)
          # Pull info from fit to find the best lambda   
          fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))
          # Needs to be min or max depending on cv measure (MSE min, AUC max, ...)
          best.lam <- fit.df[which.min(fit.df[,1]),]
          cvm.best <- best.lam[,1]
          lambda.best <- best.lam[,2]
          # Position of best lambda in cv.glmnet output
          nrow.best <- best.lam[,3]
          # Get the betas from the best lambda value
          ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
          ret[ret == 0.0] <- NA
          # Pull the non-zero betas from model
          as.vector(ret[which(!is.na(ret)),])
        },
        error = function(cond) {
          # Should fire only when all predictors have 0 variance.
          message('Error with gene ' %&% gene %&% ', index ' %&% i)
          message(geterrmessage())
          return(data.frame())
        }
      )
    } 
    model.end <- as.POSIXct(Sys.time())
    #message("        #modeling span:\t", model.end-model.start) 
    if (length(bestbetas) > 0) {   # bestbetas is a vector for non-zero betas
      #message("#end of modeling")
      #message("#number of betas:\t", length(bestbetas), '\n')
      names(bestbetas) <- rownames(ret)[which(!is.na(ret))] # names(bestbetas is a list of varID)
      # Pull out the predictions at the best lambda value.    
      pred.mat <- fit$fit.preval[,nrow.best]
      res <- summary(lm(exppheno~pred.mat))
      genename <- as.character(gene_annot[gene, 3])
      rsq <- res$r.squared
      pval <- res$coef[2,4]  
      resultsarray[gene,] <- c(gene, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval, genename, w_tf, pct_tf, ttl_snps, used_tf)  
      # Output best shrunken betas for PrediXcan
      bestbetalist <- names(bestbetas)
      bestbetainfo <- snp_annot[bestbetalist,]
      betatable <- as.matrix(cbind(bestbetainfo,bestbetas))
      betatable <- na.omit(betatable)
      #message('bestbetas=', betatable[,"varID"])
      if (nrow(betatable)>2){
      write_covariance(gene, cisgenos, betatable[,"rsid"], betatable[,"varID"], covariance_out)}
      # Output "gene", "rsid", "refAllele", "effectAllele", "beta"
      # For future: To change rsid to the chr_pos_ref_alt_build label, change "rsid" below to "varID".
      betafile <- cbind(gene,betatable[, "varID"],betatable[,"rsid"],betatable[,"refAllele"],betatable[,"effectAllele"],betatable[,"bestbetas"], alpha)
      # Transposing betafile necessary for correct output from write() function
      write(t(betafile), file = workingweight, ncolumns = 7, append = TRUE, sep = "\t")
      write(resultsarray[gene,], file = workingbest, ncolumns = 14, append = TRUE, sep = "\t")
    } else {
      #message("        #end of modeling")
      #message("        #number of betas:\t0",  '\n')
      genename <- as.character(gene_annot[gene,3])
      #resultsarray[gene,1] <- gene
      #resultsarray[gene,2:13] <- c(alpha,NA,NA,NA,0,NA,NA, genename, w_tf, pct_tf, ttl_snps, used_tf)
      #resultsarray[gene,13] <- genename
    }
  }
  write.table(resultsarray,file=out_dir %&% "TW_" %&% tis %&% "_chr" %&% chrom %&% "_exp_" %&% n_k_folds %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% ".txt",quote=F,row.names=F,sep="\t")
}

write_covariance <- function(gene, cisgenos, model_rsids, model_varIDs, covariance_out) {
  model_geno <- cisgenos[,model_varIDs, drop=FALSE]
  geno_cov <- cov(model_geno)
  cov_df <- data.frame(gene=character(),rsid1=character(),rsid2=character(), covariance=double())
  for (i in 1:length(model_rsids)) {
    for (j in i:length(model_rsids)) {
      cov_df <- rbind(cov_df, data.frame(gene=gene,rsid1=model_rsids[i], rsid2=model_rsids[j], covariance=geno_cov[i,j]))
    }
  }
  write.table(cov_df, file = covariance_out, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}
