

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
  pvals <- pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE) * 2
  
  return(list(pvals=pvals, tvals=tvals))
}


mut <- c(10,11, 12,13,15,11,12,5)

         

wt <- c(10,34, 19,18,15,22,16,78)
        


pixels <- c(mut, wt)

mat <- matrix(c(pixels, pixels), nrow = 16)


genotype <- c(rep('mutant', 8), rep('wildtype',8))

df <- data.frame(genotype=genotype)

fit <- lm(mat ~ df$genotype)

results <- pandt_vals(fit)

tvals <- results$tvals[2,]
pvals <- results$pvals[2,]

print(pvals)
print(tvals)