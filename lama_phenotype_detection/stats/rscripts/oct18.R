# 231013
# Test whether crown-rump or whole embryo volume is the best variable to account for developmental differences at E14.5
library(MASS)
library(car)


df_crl <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/staging_info_crl.csv')

df_whole <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/staging_info_volume.csv')

df_organs <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/raw_wt_organ_vols.csv')

out_dir <- '/home/neil/IMPC_research/neil/E14.5/baselines/QC/regression_plots'

individ_label_dir = file.path(out_dir, 'per_label')
dir.create(individ_label_dir, showWarnings = FALSE)

qq_dir = file.path(out_dir, 'qq')
residuals_vs_fitted_dir = file.path(out_dir, 'res_vs_fit')

dir.create(file.path(qq_dir), showWarnings = FALSE)
dir.create(residuals_vs_fitted_dir, showWarnings = FALSE)

colnames(df_organs)[1] <- 'vol'
colnames(df_whole)[1] <- 'vol'
colnames(df_whole)[2] <- 'whole'
colnames(df_crl)[2] <- 'crl'

bics_volume <- vector()
bics_boxcox <- vector()
bics_volume_crl <- vector()
bics_boxcox_crl <- vector()


r2_whole <- vector()
r2_whole_box_cox <- vector()
r2_crl <- vector()
r2_crl_boxcox <- vector()

df_test = merge(x=df_crl, y=df_whole, by.x='vol', by.y='vol')
df_test1 = merge(x=df_test, y=df_organs, by.x='vol', by.y='vol')


plot_lm <- function(fit_, title_, outpath){
  devAskNewPage(ask = FALSE)
  png(outpath, width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(fit_)
  title(title_)
  layout(1)
  dev.off()
}

for(label in names(df_test1)) {
  
  # Only organ volum columns please
  if (!startsWith(label, 'X')){
    next()
  }
  
  # No unused labels thanks
  if (all(df_test1[[label]] == 0)){
    next()
  }
  #
  # if (label != 'X187'){
  #   next()
  # }
  
  label_dir <- file.path(individ_label_dir, label)
  dir.create(label_dir, showWarnings = FALSE)
  
  print(label)
  fit_whole <- lm(df_test1[[label]] ~ df_test1$whole)
  

  Box <- boxcox(df_test1[[label]] ~ df_test1$whole, plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox = data.frame(Box$x, Box$y)
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]
  tformed <- bcPower(df_test1[[label]], lambda)

  fit_whole_box <-  lm(tformed~  df_test1$whole)
  
  r2_whole <- append(r2_whole, summary(fit_whole)$r.squared)
  r2_whole_box_cox <- append(r2_whole_box_cox, summary(fit_whole_box)$r.squared)

  bics_volume <- append(bics_volume, BIC(fit_whole))
  bics_boxcox <- append(bics_boxcox, BIC(fit_whole_box))


  ################################## whole embryo
  ###### no normalisation

  ## Save qq plots
  title_ = paste(label, 'whole embryo')
  plot_path = file.path(label_dir, paste(label, '.png'))
  plot_lm(fit_whole, title_, plot_path)

  #####  boxcox transformed
  ## Save qq plots
  # devAskNewPage(ask = FALSE)
  title_ = paste(label, 'whole embryo boxcox')
  plot_path = file.path(label_dir, paste(label, 'boxcox.png'))
  plot_lm(fit_whole_box, title_, plot_path)

  ############## Using CRl ################################
  
  fit_crl <- lm(log(df_test1[[label]]) ~ df_test1$crl)
  
  Box <- boxcox(df_test1[[label]] ~ df_test1$crl, plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox <- data.frame(Box$x, Box$y)
  Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),]
  lambda <- Cox2[1, "Box.x"]
  tformed <- bcPower(df_test1[[label]], lambda)
  
  fit_crl_box <-  lm(tformed~  df_test1$crl)

  r2_crl <- append(r2_crl, summary(fit_crl)$r.squared)
  r2_fit_crl_box <- append(r2_crl_boxcox, summary(fit_crl_box)$r.squared)
  
  bics_volume_crl <- append(bics_volume_crl, BIC(fit_crl))
  bics_boxcox_crl <- append(bics_boxcox_crl, BIC(fit_crl_box))

  ## Save qq plots
  plot_path = file.path(label_dir, paste(label, 'crl.png'))
  title_ <- paste(label, 'CRL')
  plot_lm(fit_crl, title_, plot_path)
  
  #####  boxcox transformed
  ## Save qq plots
  plot_path = file.path(label_dir, paste(label, 'boxcox_crl.png'))
  title_ <- paste(label, 'CRL boxcox')
  plot_lm(fit_crl_box, title_, plot_path)


}



## Save bics plots etc
bics_df <- data.frame(bics_volume, bics_boxcox, bics_volume_crl, bics_boxcox_crl)
colnames(bics_df) <- c('whole_embryo', 'whole_embryo_bc', 'crl', 'crl_bc')

write.csv(bics_df, file.path(out_dir, 'bics.csv'))

boxplot_path = file.path(out_dir, 'bics.png')
boxplot(bics_df)
png(file=boxplot_path)
dev.off()


## save r2 plots and csv
r2_df <- data.frame(r2_whole, r2_whole_box_cox, r2_crl, r2_crl_boxcox)
colnames(bics_df) <- c('whole_embryo', 'whole_embryo_bc', 'crl', 'crl_bc')

write.csv(r2_df, file.path(out_dir, 'r2.csv'))

boxplot_path = file.path(out_dir, 'r2.png')
boxplot(bics_df)
png(file=boxplot_path)
dev.off()




