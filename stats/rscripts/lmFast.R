library(MASS)
library(car)

args <- commandArgs(trailingOnly = TRUE);


testing = FALSE;

if (testing == FALSE){
  pixels_file <- args[1];  # A binary containing the voxel to be tested. Masked voxels will have been removed
  groups_file <- args[2];  # CSV containing the genotype and crown-rum (or other staging metric)
  pvals_out <- args[3];    # The output file path for the pvalues
  tvals_out <- args[4];    # The output file path for the t-statistics
  formula <- args[5];      # The formula to use.
  do_box_cox <- args[6];      # The formula to use.
}else{
  pixels_file <- 'test_data_for_R_LM/testpixelfile';
  groups_file <- 'test_data_for_R_LM/groups.csv';
  pvals_out <- "~/test_pvals.bin";
  tvals_out <- "~/test_tscores.bin";
  formula <- "genotype";
  do_box_cox <- TRUE;
}

# Create a data frame of the groups
g <- read.table(groups_file, header=TRUE, sep=',')
groups <- data.frame(g)

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
  print(typeof(tvals))
  pvals <- pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE) * 2
  
  return(list(pvals=pvals, tvals=tvals))
}

boxy <- function(single_organ_data, row_indices){
  # Do a boxcox tranformon the data
  # If row_indices subset based on threse rows (when doing specimen n =1)
 
  if (identical(row_indices, FALSE)){
    Box <- boxcox(single_organ_data ~ groups$crl, plotit = FALSE, lambda = seq(-2, 2, len = 1000))
  }else{
    single_organ_data <- single_organ_data[row_indices]
    Box <- boxcox(single_organ_data ~ groups$crl[row_indices], plotit = FALSE, lambda = seq(-2, 2, len = 1000))
  }

  Cox = data.frame(Box$x, Box$y)  
  CoxSorted = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = CoxSorted[1, "Box.x"]    
  tformed <- bcPower(single_organ_data, lambda)
  return(tformed)
}

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- abs(matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2]))
close(con)


formula_elements <- strsplit(formula, split=',')
print('lm formula elements');
print(formula_elements)



if (do_box_cox == TRUE){
  print('##doing boxcox##')
  tformed = apply(mat, 2, boxy, row_indices=FALSE)
  fit <- lm(tformed ~., data=groups[, unlist(formula_elements)])

}else{
  fit <- lm(mat ~., data=groups[, unlist(formula_elements)])
}

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
  
  spec_results <- pandt_vals(fit_specimen)
  pvals = append(pvals, spec_results$pvals[2,])
  tscores = append(tscores, spec_results$tvals[2,])

}

print(length(pvals))
poutCon <- file(pvals_out, "wb")
# writeBin(results$pvals[2,], poutCon)
writeBin(pvals, poutCon)
close(poutCon)

toutCon <- file(tvals_out, "wb")
# writeBin(0 - results$tvals[2,], toutCon)
writeBin(0 - tscores, toutCon)


close(toutCon)


