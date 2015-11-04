args <- commandArgs(trailingOnly = TRUE);
pixels_file <- args[1];
groups_file <- args[2];
pvals_out <- args[3];
tvals_out <- args[4];
#pixels_file <- '/tmp/raw_data_for_r.csv';
#groups_file <- '/tmp/groups_for_liear_model.csv';

pixels <- read.csv(pixels_file, sep=',',header=FALSE)
groups <- read.csv(groups_file, header=TRUE)




df <- data.frame(pixels, groups)

pvals <- c()
tvals <- c()

for (i in 2:length(df) -1){
  lm1 = lm(df[,i] ~ df$sex);
  pval <- summary(lm1)$coef[,"Pr(>|t|)"][2];
  tval <- summary(lm1)$coef[,"t value"][2];
  pvals[i] <- pval;
  tvals[i] <- tval
};

write.csv(pvals, file=pvals_out, row.names=FALSE, col.names=FALSE);
write.csv(tvals, file=tvals_out, row.names=FALSE, col.names=FALSE);

