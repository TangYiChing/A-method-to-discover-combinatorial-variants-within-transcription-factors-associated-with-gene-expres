

start_time <- as.POSIXct(Sys.time())
TW_CV_model(expression_RDS, genotype_txt, gene_annot_RDS, snp_annot_RDS, tf_gene_annot_RDS,
    n_k_folds, alpha, out_dir, study, chrom, snpset, window, tf_genotype_dir, tf_snp_annot_dir, seed=seed)
end_time <- as.POSIXct(Sys.time())
time.taken <- end_time - start_time
message("Program Start: \t", start_time)
message("Program End: \t", end_time)
print(time.taken)
