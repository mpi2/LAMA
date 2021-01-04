# TODO: The specimen-level ananlysis is only done if there if the groups file has a genotype column
# Containing 'mutant' and 'wildtype' entries. For the permutation stats, these do not exist
# This is fine as we don't want specimen-level calls for permutation testing. Make this more logical.

library(MASS)


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

# boxy <- function(single_organ_data, row_indices){
#   # Do a boxcox tranformon the data
#   # If row_indices subset based on these rows (when doing specimen n =1)
#
#   if (identical(row_indices, FALSE)){
#     Box <- boxcox(single_organ_data ~ groups$crl, plotit = FALSE, lambda = seq(-2, 2, len = 1000))
#   }else{
#     single_organ_data <- single_organ_data[row_indices]
#     Box <- boxcox(single_organ_data ~ groups$crl[row_indices], plotit = FALSE, lambda = seq(-2, 2, len = 1000))
#   }
#
#   Cox = data.frame(Box$x, Box$y)
#   CoxSorted = Cox[with(Cox, order(-Cox$Box.y)),]
#   lambda = CoxSorted[1, "Box.x"]
#   tformed <- bcPower(single_organ_data, lambda)
#   return(tformed)
# }

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- abs(matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2]))
close(con)


formula_elements <- strsplit(formula, split=',')
# print('lm formula elements');
# print(formula_elements);



if (do_box_cox == TRUE){
  print('##doing boxcox##')
  # tformed = apply(mat, 2, boxy, row_indices=FALSE)
  # fit <- lm(tformed ~., data=groups[, unlist(formula_elements)])

}else{
  fit <- lm(mat ~., data=groups[, unlist(formula_elements)])
}

# line_level_plot_dir <-  file.path(plot_dir, 'line_level_plots')
# dir.create(line_level_plot_dir, showWarnings = FALSE)

# apply(mat, 2, plot_lm, groups, outdir=line_level_plot_dir)

results <- pandt_vals(fit)
pvals = results$pvals[2,]
tscores = results$tvals[2,]


# Now fit each specimen individually to the linear model
mutant_row_nums = which(groups$genotype == 'mutant');
wt_row_nums = which(groups$genotype == 'wildtype')

for (r in mutant_row_nums){
  #For each mutant add the mutant row number to the wt row indices
  row_indices = c(wt_row_nums, r)

  if (do_box_cox == TRUE){

    tformed = apply(mat, 2, boxy, row_indices=row_indices)
    fit_specimen <- lm(tformed ~., data=groups[row_indices, unlist(formula_elements)])

  }else{
    fit_specimen <- lm(mat[row_indices, ] ~., data=groups[row_indices, unlist(formula_elements)])

  }

  specimen_plot_path = file.path(plot_dir, groups[r, 0])
  #plot_lm(specimen_plot_path, fit_specimen)

  spec_results <- pandt_vals(fit_specimen)
  pvals = append(pvals, spec_results$pvals[2,])
  tscores = append(tscores, spec_results$tvals[2,])

}

# print('writigpvals file')
poutCon <- file(pvals_out, "wb")
# writeBin(results$pvals[2,], poutCon)
writeBin(pvals, poutCon)
close(poutCon)

toutCon <- file(tvals_out, "wb")

# R returns the genotype effect for wildtype so we must flip the sign to get it for mutant
writeBin(0 - tscores, toutCon)


close(toutCon)


