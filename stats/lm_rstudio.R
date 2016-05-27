

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


mut <- c(-0.07058107, -0.05060999, -0.05551733, -0.04817696, -0.05799195)
wt <- c(-0.04575749, -0.02227639, -0.02687215, -0.01322827, -0.01772877)


#wt <- c(rep(37, 11))
#mut <- c(rep(38, 3))

pixels <- c(mut, wt)

mat <- matrix(c(pixels, pixels), nrow = 10)


genotype <- c(rep('mutant', 5), rep('wildtype',5))
#sex <- c('m', 'f', 'm', 'f', 'm', 'f', 'm', 'm', 'm', 'f', 'm', 'f', 'f', 'f')
#sex <- c(rep('m', 11), rep('f', 3))

df <- data.frame(genotype=genotype)

fit <- lm(mat ~ df$genotype)

results <- pandt_vals(fit)

tvals <- results$tvals[2,]
pvals <- results$pvals[2,]

summary(fit)