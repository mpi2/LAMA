#args <- commandArgs(trailingOnly = TRUE);

#organ_vol_file <- args[1];
#groups_file <- args[2];
#out_file <- args[3];

organ_vol_file <- '/home/neil/share/test_dataset_for_lama_dev/mut/output/stats/organ_vols_combined.csv'
groups_file <- '/home/neil/share/test_dataset_for_lama_dev/mut/output/stats/groups_for_vol.csv'
out_file <- '/home/neil/share/test_dataset_for_lama_dev/mut/output/stats/organ_vol_stats'
formula <- 'genotype';
formula_elements <- strsplit(formula, split=',')

vols <- as.matrix(read.table(organ_vol_file, sep=',', header=FALSE))
groups <- data.frame(read.table(groups_file, sep=',', header=TRUE))

fit <- lm(vols ~., data=groups[, unlist(formula_elements)])
summary(fit)