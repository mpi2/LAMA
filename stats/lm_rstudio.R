
pixels_file <- '/tmp/tmp_data_for_lm';
groups_file <- '/home/neil/share/test_dataset_for_lama_dev/mut/output/stats/combined_groups.csv'
pvals_out <- '/tmp/pvals_out.dat'
formula <- 'genotype'; # Just foing one formula at a time for now

# Create a data frame of the groups
g <- read.table(groups_file, header=TRUE, sep=',')
groups <- data.frame(g)



pvalOnly2 <- function(fit) {
  # from: http://stackoverflow.com/questions/33652502/quickly-retrieve-pvalues-from-multiple-lm-in-r/33664809#33664809
  # get estimates
  # This gets the pvalues for each fixed effect in the model
  est <- fit$coefficients[fit$qr$pivot, ]
  
  # get R: see stats:::summary.lm to see how this is calculated
  p1 <- 1L:(fit$rank)
  R <- diag(chol2inv(fit$qr$qr[p1, p1, drop = FALSE]))
  
  # get residual sum of squares for each
  resvar <- colSums(fit$residuals^2) / fit$df.residual
  # R is same for each coefficient, resvar is same within each model 
  se <- sqrt(outer(R, resvar))
  
  pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE) * 2
}

modelPvalOnly <- function(fit) {
  # This get the p-value for the whole model
  f <- t(fit$fitted.values)
  if (attr(fit$terms, "intercept"))  {
    mss <- rowSums((f - rowMeans(f)) ^ 2)
    numdf <- fit$rank - 1
  } else {
    mss <- rowSums(f ^ 2)
    numdf <- fit$rank
  }
  
  resvar <- colSums(fit$residuals^2) / fit$df.residual
  fstat <- (mss / numdf) / resvar
  pval <- pf(fstat, numdf, fit$df.residual, lower.tail = FALSE)
  pval
}

con <- file(pixels_file, "rb")
dim <- readBin(con, "integer", 2)
mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)


formula_elements <- strsplit(formula, split=' ')
fit <- lm(mat ~., data=groups[, unlist(formula_elements)])

p <- modelPvalOnly(fit)

#toutCon <- file(tvals_out, "wb")
#writeBin(tvals, toutCon)
#close(toutCon)

poutCon <- file(pvals_out, "wb")
writeBin(p, poutCon)
close(poutCon)



