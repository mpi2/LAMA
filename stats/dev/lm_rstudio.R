

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


mut <- c(1,
         2,
         3,
         4
         
)
wt <- c(5,
        6,
        7,
        8
        
)

pixels <- c(mut, wt)

mat <- matrix(c(pixels, pixels), nrow = 8)


genotype <- c(rep('mutant', 4), rep('wildtype',4))
#crl <- c(6,7,9,7, 2,3,2,3)

df <- data.frame(genotype=genotype)

#fit_with_crl <- lm(mat ~ df$genotype+df$crl)
fit_no_crl <- lm(mat ~ df$genotype)


results <- pandt_vals(fit_no_crl)

tvals <- results$tvals[2,]
pvals <- results$pvals[2,]

print(pvals)
print(tvals)
