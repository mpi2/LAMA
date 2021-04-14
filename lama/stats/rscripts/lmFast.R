# TODO: The specimen-level ananlysis is only done if there if the groups file has a genotype column
# Containing 'mutant' and 'wildtype' entries. For the permutation stats, these do not exist
# This is fine as we don't want specimen-level calls for permutation testing. Make this more logical.


library(MASS)

#install.packages('abind', repos='http://cran.us.r-project.org')
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

 
aov_fstat_pvals <- function(fit) {
    
    # str(fit) - #uncomment this line to see what fit$assign etc is if you're having trouble with understanding the code
    
    # yes I hate this too - but this is the only way to get the proper p-values
    # I had to bascially integrate code from aov.summary to get the interaction p-vals     
    
    # get the index values for each set of p-values 
    asgn <- fit$assign[fit$qr$pivot[1L:fit$rank]]
     
    # this just gets the columns of interest - remove the first value as its just the intercept - which we don't care about
    # additionally remove the fourth column, which should be staging. (the genotype:treatment value is added last in the aov, irregardless of the call...)
    if (length(asgn) == 5) { 
	    asgn <- asgn[-4]
	    asgn <-asgn[-1]
    
    } else { # these are one-ways, which should just have either the genotype/treatment effect in the second column
	    assgn <- asgn[2]
    }

    uasgn <- unique(asgn)
    nterms <- length(uasgn)
    effects <- fit$effects 
  
    if(!is.null(effects)) {
	effects <- as.matrix(effects)[seq_along(asgn),,drop=FALSE]      
    }
    # remove undesired columns same as asgn
    if (NCOL(effects) == 5) {
        effects <- effects[,-4]
        effects <- effects[,-1]
    }
    rdf <- fit$df.residual

    resid <- as.matrix(fit$residuals)
    
    nresp <- NCOL(resid)  
             
    df <- ss <- matrix(, nrow = nresp, ncol = (nterms+1)) # need extra column for residuals
  
    for (y in 1L:nresp) {
        
        for(i in seq(nterms)) { 		
            ai <- (asgn == uasgn[i])
            
            df[y, i] <- sum(ai) # this may be botch as I'm removing asgn columns

            ss[y, i] <- sum(effects[ai, y]^2)
	        
	}

    }
    

    if(rdf > 0L) {
	df[, (nterms+1)] <- rdf
        ss[, (nterms+1)] <- sum(resid[, y]^2)
	#nmrows <- c(nmrows,  "Residuals")
     }
    
    nt <- NCOL(df)
    ms <- ifelse(df > 0L, ss/df, NA)
    
    
    if(rdf > 0L) {
	fstat <- ms/ms[nt] # divide all mean squares by the residuals (last column)
        
        fstat <- fstat[,-nt] # we don't want the last column anymore 
	
                
         
        pvals <- pf(fstat, df, rdf, lower.tail = FALSE) 
	# This adds the residuals so you need to remove them again?
        pvals <- pvals[,-nt]
	} 
    return(list(pvals=pvals, tvals=fstat)) 

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
    fit <- aov(mat ~genotype:treatment+., data=groups[, unlist(formula_elements)])
    results <- aov_fstat_pvals(fit)
    dim <- c(length(results$pval[,1]), length(results$pval[1,]), 0)
    fscores <- pvals <- array(numeric(), dim)
    pvals <- abind(pvals, data.matrix(results$pvals))
    # for g by e studies, the aov returns f_values
    fscores <- abind(fscores, data.matrix(results$tvals))
    
    

} else {
    fit <- lm(mat ~., data=groups[, unlist(formula_elements)])
    results <- pandt_vals(fit)
    pvals = results$pvals[2,]
    tscores = results$tvals[2,]} 
       
# line_level_plot_dir <-  file.path(plot_dir, 'line_level_plots')
# dir.create(line_level_plot_dir, showWarnings = FALSE)

# apply(mat, 2, plot_lm, groups, outdir=line_level_plot_dir)





# Now fit each specimen individually to the linear model
if("treatment" %in% colnames(groups)){
    mutant_row_nums = which((groups$genotype == 'mutant') | (groups$treatment == 'treatment'))
    wt_row_nums = which((groups$genotype == 'wildtype') & (groups$treatment == 'vehicle'))
    non_interaction_row_nums = which(!((groups$genotype == 'mutant') & (groups$treatment == 'treatment')))

}else{
    mutant_row_nums = which(groups$genotype == 'mutant');
    wt_row_nums = which(groups$genotype == 'wildtype')
}
 
for (r in mutant_row_nums){
  #For each mutant add the mutant row number to the wt row indices
  print(groups[r,])
  row_indices = c(wt_row_nums, r)
  # I'm just going to comment this out to improve runtime
  # if (do_box_cox == TRUE){

  # tformed = apply(mat, 2, boxy, row_indices=row_indices)
  #  fit_specimen <- lm(tformed ~., data=groups[row_indices, unlist(formula_elements)])

  # if ("treatment" %in% colnames(groups)){
  #  row_elements = c(unlist(formula_elements))
    # if the row contain's wild-type, it means you should only test treatment
    # if (groups[r,]$genotype == "wildtype") {
    #	row_elements <- row_elements[! row_elements %in% c("genotype")]
    #    fit_specimen <- aov(mat[row_indices, ] ~., data=groups[row_indices, unlist(row_elements)])

    # spec_results <- aov_fstat_pvals(fit_specimen)
        # create voxel_no X spec_num X 3 array of NULLs
     #   dim <- c(length(spec_results$pval[,1]), length(results$pval[1,]), 3)
     # spec_fstat <- spec_pvals <- array(data = NA, dim)
        
	# spec_fstat[ , ,2] = data.matrix(spec_results$tvals)

	# spec_pvals[ , ,2] = data.matrix(spec_results$pvals)
        # fscores = abind(fscores, spec_fstat)
	# pvals = abind(pvals, spec_pvals)
        
  

        #as this is treatment only, the second column in the dataset should contian the t-stats and pvals 
    #}
    # if the row contains vehicle, it means you should only test genotype
    #if (groups[r,]$treatment == "vehicle") {
    #    row_elements <- row_elements[! row_elements %in% c("treatment")]
    #    fit_specimen <- aov(mat[row_indices, ] ~., data=groups[row_indices, unlist(row_elements)])
        
    #	spec_results <- aov_fstat_pvals(fit_specimen)

	
	# create voxel_no X spec_num X 3 array of NULLs
     #  dim <- c(length(spec_results$pval[,1]), length(results$pval[1,]), 3)
	# spec_fstat <- spec_pvals <- array(data = NA, dim)
        
	# spec_fstat[ , ,1] = data.matrix(spec_results$tvals)
        # spec_pvals[ , ,1] = data.matrix(spec_results$pvals)
	
	# fscores = abind(fscores, spec_fstat)
        # pvals = abind(pvals, spec_pvals)

    #} 
    # these are the g-by-e interaction specimens, perform a two-way anova on them
    if ((groups[r,]$genotype == 'mutant') & (groups[r,]$treatment == 'treatment')) {
        row_indices = c(non_interaction_row_nums, r)
        fit_specimen <- aov(mat[row_indices, ] ~genotype:treatment+ ., data=groups[row_indices, unlist(row_elements)])
	
	spec_results <- aov_fstat_pvals(fit_specimen)

	pvals = abind(pvals, data.matrix(spec_results$pvals))
	fscores = abind(fscores, data.matrix(spec_results$tvals))
	
    }
    
    
  } else {
    fit_specimen <- lm(mat[row_indices, ] ~., data=groups[row_indices, unlist(formula_elements)])

    spec_results <- pandt_vals(fit_specimen)
    pvals = append(pvals, spec_results$pvals[2,])
    tscores = append(tscores, spec_results$tvals[2,])
  }
  
  specimen_plot_path = file.path(plot_dir, groups[r, 0])
  #plot_lm(specimen_plot_path, fit_specimen)


        
}
# print('writigpvals file')

if ("treatment" %in% colnames(groups)) {
    file_exts = c('genotype','treatment','interaction')
    for (f in 1:length(results$pval[1,])) {
	
        f_data <- fscores[ ,f, ]  
        p_data <- pvals[ ,f, ]
                
        pvals_g_out = paste(pvals_out, file_exts[f], sep="_")
	
	poutCon <- file(pvals_g_out, "wb")
	    
	writeBin(as.vector(p_data), poutCon)
	close(poutCon)
        
	fvals_g_out = paste(tvals_out, file_exts[f], sep="_")
	foutCon <- file(fvals_g_out, "wb")

        # F stat measues if the means between groups are equal or not, it has no direction. 
	# Therefore you don't need to flip the sign
	writeBin(as.vector(f_data), foutCon)
        close(foutCon)
    }
        
   
} else {

    poutCon <- file(pvals_out, "wb")
    # writeBin(results$pvals[2,], poutCon)
    writeBin(pvals, poutCon)
    close(poutCon)

    toutCon <- file(tvals_out, "wb")

    # R returns the genotype effect for wildtype so we must flip the sign to get it for mutant
    writeBin(0 - tscores, toutCon)


    close(toutCon)

    }
