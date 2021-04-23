# TODO: The specimen-level ananlysis is only done if there if the groups file has a genotype column
# Containing 'mutant' and 'wildtype' entries. For the permutation stats, these do not exist
# This is fine as we don't want specimen-level calls for permutation testing. Make this more logical.

library(MASS)

# install the abind package if not previously found

if (!require(abind)) install.packages('abind', repos='http://cran.us.r-project.org')
library(abind)

args <- commandArgs(trailingOnly = TRUE);

pixels_file <- args[1];  # A binary containing the voxel to be tested. Masked voxels will have been removed
groups_file <- args[2];  # CSV containing the genotype and crown-rum (or other staging metric)
pvals_out <- args[3];    # The output file path for the pvalues
tvals_out <- args[4];    # The output file path for the t-statistics
formula <- args[5];      # The formula to use.
do_box_cox <- args[6];      # The formula to use.
plot_dir <- args[7];

# Create a data frame of the groups
g <- read.table(groups_file, header=TRUE, sep=',')
groups <- data.frame(g)

counter = 0

plot_lm <- function(data, groups, outdir){
  return()
  counter <- counter + 1;
  outname = file.path(outdir, paste(counter, '.png'))
  png(outname, width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  fit_plot <- lm(data ~ groups$crl)
  plot(fit_plot)
#  bic <- BIC(fit_plot)
#  print(bic)
  layout(1)
  dev.off()
}

pandt_vals <- function(fit) {
  # get estimates
  est <- fit$coefficients[fit$qr$pivot, ]

  # get R: see stats:::summary.lm to see how this is calculated
  p1 <- 1L:(fit$rank)
  R <- diag(chol2inv(fit$qr$qr[p1, p1, drop = FALSE]))

  # get residual sum of squares for each
  resvar <- colSums(fit$residuals^2) / fit$df.residual
  # R is same for each coefficient, resvar is same within each model
  se <- sqrt(outer(R, resvar))

  tvals <- est / se
  #print(typeof(tvals))
  pvals <- pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE) * 2

  return(list(pvals=pvals, tvals=tvals))
}

con <- file(pixels_file, "rb")

dim <- readBin(con, "integer", 2)
mat <- abs(matrix(readBin(con, "numeric", prod(dim)), dim[1], dim[2]))

close(con)


formula_elements <- strsplit(formula, split=',')

# just commenting out to improve speed
# if (do_box_cox == TRUE){
  # print('##doing boxcox##')
  # tformed = apply(mat, 2, boxy, row_indices=FALSE)
  # fit <- lm(tformed ~., data=groups[, unlist(formula_elements)])

if ("treatment" %in% colnames(groups)){
    fit <- lm(mat ~genotype:treatment+., data=groups[, unlist(formula_elements)])

    results <- pandt_vals(fit)

    tval <- results$tvals[c(2,3,5),]
    pval <- results$pvals[c(2,3,5),]

    dim <- c(length(pval[,1]), length(pval[1,]), 0)


    tscores <- pvals <- array(numeric(), dim)

    pvals <- abind(pvals, data.matrix(pval))
    # for g by e studies, the aov returns f_values
    tscores <- abind(tscores, data.matrix(tval))



} else {
    fit <- lm(mat ~., data=groups[, unlist(formula_elements)])
    results <- pandt_vals(fit)
    pvals = results$pvals[1,]
    tscores = results$tvals[2,]}



# Now fit each specimen individually to the linear model
if("treatment" %in% colnames(groups)){
    # I just want the interaction rows to reduce multiple testing

    mutant_row_nums = which((groups$genotype == 'mutant') & (groups$treatment == 'treatment'))
    # wt_row_nums = which((groups$genotype == 'wildtype') & (groups$treatment == 'vehicle'))
    non_interaction_row_nums = which(!((groups$genotype == 'mutant') & (groups$treatment == 'treatment')))

    for (r in mutant_row_nums){
        #For each mutant add the mutant row number to the wt row indices
        row_indices = c(non_interaction_row_nums, r)

        fit_specimen <- lm(mat[row_indices, ] ~genotype:treatment+ ., data=groups[row_indices, unlist(formula_elements)])

        spec_results <- pandt_vals(fit_specimen)

        pval <- results$pvals[c(2,3,5),]

        tval <- results$tvals[c(2,3,5),]

        pvals = abind(pvals, data.matrix(pval))

        tscores = abind(tscores, data.matrix(tval))
	    
        specimen_plot_path = file.path(plot_dir, groups[r, 0])
    }
    #TODO: probably don't hard code file extensions
    file_exts = c('genotype','treatment','interaction')
    for (f in 1:3) {
	    
        t_data <- tscores[f, , ]
        p_data <- pvals[f, , ]

        pvals_g_out = paste(pvals_out, file_exts[f], sep="_")


       
      
        poutCon <- file(pvals_g_out, "wb")
        writeBin(as.vector(p_data), poutCon)
        close(poutCon)

        tvals_g_out = paste(tvals_out, file_exts[f], sep="_")
        toutCon <- file(tvals_g_out, "wb")

        writeBin(0 - as.vector(t_data), toutCon)
        close(toutCon)
    }

}else{
    mutant_row_nums = which(groups$genotype == 'mutant');
    wt_row_nums = which(groups$genotype == 'wildtype')

    for (r in mutant_row_nums){

        row_indices = c(wt_row_nums, r)

        fit_specimen <- lm(mat[row_indices, ] ~., data=groups[row_indices, unlist(formula_elements)])

        spec_results <- pandt_vals(fit_specimen)
        pvals = append(pvals, spec_results$pvals[2,])
        tscores = append(tscores, spec_results$tvals[2,])

        specimen_plot_path = file.path(plot_dir, groups[r, 0])
    }
    poutCon <- file(pvals_out, "wb")
    # writeBin(results$pvals[2,], poutCon)
    writeBin(pvals, poutCon)
    close(poutCon)

    toutCon <- file(tvals_out, "wb")

    # R returns the genotype effect for wildtype so we must flip the sign to get it for mutant
    writeBin(0 - tscores, toutCon)

    close(toutCon)


}
