# 231013
# Test whether crown-rump or whole embryo volume is the best variable to account for developmental differences at E14.5
library(MASS)
library(car)


df_crl <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/staging_info_crl.csv')

df_whole <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/staging_info_volume.csv')

df_organs <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/raw_wt_organ_vols.csv')

out_dir <- '/home/neil/IMPC_research/neil/E14.5/baselines/QC/regression_plots'

individ_label_dir = file.path(out_dir, 'per_label')

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


df_test = merge(x=df_crl, y=df_whole, by.x='vol', by.y='vol')
df_test1 = merge(x=df_test, y=df_organs, by.x='vol', by.y='vol')

for(label in names(df_test1)) {
  
  # Only organ volum columns please
  if (!startsWith(label, 'X')){
    next()
  }
  
  # No unused labels thanks
  if (all(df_test1[[label]] == 0)){
    next()
  }
  
  if (label != 'X27'){
    next()
  }
  
  label_dir <- file.path(individ_label_dir, label)
  dir.create(label_dir, showWarnings = FALSE)
  
  print(label)
  fit_whole <- lm(log(df_test1[[label]]) ~ log(df_test1$whole))
  

  Box <- boxcox(fit_whole, plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox = data.frame(Box$x, Box$y)  
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]    
  tformed <- bcPower(df_test1[[label]], lambda)
  
  fit_whole_box <-  lm(tformed~  log(df_test1$whole))
  
  bics_volume <- append(bics_volume, BIC(fit_whole))
  bics_boxcox <- append(bics_boxcox, BIC(fit_whole_box))
  
  
  ################################## whole embryo
  ###### no normalisation
  
  ## Save qq plots
  #devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, '.png'))
  #png(file=qq_plot_path)
  plot(fit_whole)
  title(paste(label, 'whole embryo'))
  #dev.off()
  
  
  # resifdual vs fitted
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'res_vs_fit.png'))
  png(file=qq_plot_path)
  plot(fit_whole, 1)
  title(paste(label, 'whole embryo'))
  dev.off()
  

  # scale/location
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'scale_location.png'))
  png(file=qq_plot_path)
  plot(fit_whole, 3)
  title(paste(label, 'whole embryo'))
  dev.off()
  
  
  
  #####  boxcox transformed
  ## Save qq plots
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'boxcox.png'))
  png(file=qq_plot_path)
  plot(fit_whole_box, 2)
  title(paste(label, 'whole embryo boxcox'))
  dev.off()
  
  
  # resisdual vs fitted
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'res_vs_fit_boxcox.png'))
  png(file=qq_plot_path)
  plot(fit_whole_box, 1)
  title(paste(label, 'whole embryo boxcox'))
  dev.off()
  
  # scale/location
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'scale_location_box_cox.png'))
  png(file=qq_plot_path)
  plot(fit_whole_box, 3)
  title(paste(label, 'whole embryo boxcox'))
  dev.off()
  
  
  ############## Using CRl ################################
  
  fit_crl <- lm(log(df_test1[[label]]) ~ df_test1$crl)
  
  Box <- boxcox(log(df_test1[[label]]) ~ df_test1$crl, plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox = data.frame(Box$x, Box$y)  
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]    
  tformed <- bcPower(df_test1[[label]], lambda)
  
  fit_crl_box <-  lm(tformed~  df_test1$crl)
  
  bics_volume_crl <- append(bics_volume_crl, BIC(fit_crl))
  bics_boxcox_crl <- append(bics_boxcox_crl, BIC(fit_crl_box))
  
  
  
  
  
  ## Save qq plots
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'crl.png'))
  png(file=qq_plot_path)
  plot(fit_crl, 2)
  title(paste(label, 'CRL'))
  dev.off()
  
  
  # resifdual vs fitted
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'res_vs_fit_crl.png'))
  png(file=qq_plot_path)
  plot(fit_crl, 1)
  title(paste(label, 'CRL'))
  dev.off()
  
  
  # scale/location
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'scale_location_crl.png'))
  png(file=qq_plot_path)
  plot(fit_crl, 3)
  title(paste(label, 'CRL'))
  dev.off()
  
  
  
  #####  boxcox transformed
  ## Save qq plots
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'boxcox_crl.png'))
  png(file=qq_plot_path)
  plot(fit_crl_box, 2)
  title(paste(label, 'CRL boxcox'))
  dev.off()
  
  
  # resisdual vs fitted
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'res_vs_fit_boxcox_crl.png'))
  png(file=qq_plot_path)
  plot(fit_crl_box, 1)
  title(paste(label, 'CRL boxcox'))
  dev.off()
  
  # scale/location
  devAskNewPage(ask = FALSE)
  qq_plot_path = file.path(label_dir, paste(label, 'scale_location_box_cox_crl.png'))
  png(file=qq_plot_path)
  plot(fit_crl_box, 3)
  title(paste(label, 'CRL boxcox'))
  dev.off()
  
  
}

bics_df <- data.frame(bics_volume, bics_boxcox, bics_volume_crl, bics_boxcox_crl)
colnames(bics_df) <- c('nonorm', 'boxcox', 'nonorm_crl', 'boxcox_crl')

write.csv(bics_df, file.path(out_dir, 'bics.csv'))

boxplot_path = file.path(out_dir, 'bics.png')
boxplot(bics_df)
png(file=boxplot_path)
dev.off()





