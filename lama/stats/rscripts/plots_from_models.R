# 231013
# Test whether crown-rump or whole embryo volume is the best variable to account for developmental differences at E14.5
library(MASS)
library(car)


df_crl <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/staging_info_crl.csv')

df_whole <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/staging_info_volume.csv')

df_organs <- read.csv('/home/neil/IMPC_research/neil/E14.5/baselines/raw_wt_organ_vols.csv')

out_dir <- '/home/neil/IMPC_research/neil/E14.5/baselines/QC/regression_plots_311018'
dir.create(out_dir, showWarnings = FALSE)

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




r2_whole <- vector()
r2_whole_log <- vector()
r2_whole_box_cox <- vector()
r2_whole_log_cox <- vector()
r2_crl <- vector()
r2_crl_boxcox <- vector()
labels <- vector()

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

  label_dir <- file.path(individ_label_dir, label)
  dir.create(label_dir, showWarnings = FALSE)
  
  print(label)
  labels <- append(labels, label)



  ################################## whole embryo########################################
  fit_whole <- lm(df_test1[[label]] ~ df_test1$whole)
  fit_whole_log <- lm(log(df_test1[[label]]) ~ log(df_test1$whole))
  

  Box <- boxcox(df_test1[[label]] ~ df_test1$whole, plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox = data.frame(Box$x, Box$y)
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]
  tformed <- bcPower(df_test1[[label]], lambda)

  Box <- boxcox(log(df_test1[[label]]) ~ log(df_test1$whole), plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox = data.frame(Box$x, Box$y)
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
  lambda = Cox2[1, "Box.x"]
  tformed_log <- bcPower(df_test1[[label]], lambda)

  fit_whole_box <-  lm(tformed~  df_test1$whole)
  fit_whole_log_box <-  lm(tformed_log~  df_test1$whole)
  
  r2_whole <- append(r2_whole, summary(fit_whole)$r.squared)
  r2_whole_log <- append(r2_whole_log, summary(fit_whole_log)$r.squared)
  r2_whole_box_cox <- append(r2_whole_box_cox, summary(fit_whole_box)$r.squared)
  r2_whole_log_cox <- append(r2_whole_log_cox, summary(fit_whole_log_box)$r.squared)


  ## Save qq plots
  title_ = paste(label, 'whole embryo')
  plot_path = file.path(label_dir, paste(label, 'whole.png'))
  plot_lm(fit_whole, title_, plot_path)

  title_ = paste(label, 'whole embryo boxcox')
  plot_path = file.path(label_dir, paste(label, 'whole_boxcox.png'))
  plot_lm(fit_whole_box, title_, plot_path)

  title_ = paste(label, 'whole embryo log_boxcox')
  plot_path = file.path(label_dir, paste(label, 'whole_boxcox_log.png'))
  plot_lm(fit_whole_log_box, title_, plot_path)

  title_ = paste(label, 'whole embryo log')
  plot_path = file.path(label_dir, paste(label, 'whole_log.png'))
  plot_lm(fit_whole_log, title_, plot_path)




  ############## Using CRl ################################
  
  fit_crl <- lm(log(df_test1[[label]]) ~ df_test1$crl)

  Box <- boxcox(df_test1[[label]] ~ df_test1$crl, plotit = FALSE, lambda = seq(-2, 2, len = 4000))
  Cox <- data.frame(Box$x, Box$y)
  Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),]
  lambda <- Cox2[1, "Box.x"]
  tformed <- bcPower(df_test1[[label]], lambda)

  fit_crl_box <-  lm(tformed~  df_test1$crl)

  r2_crl <- append(r2_crl, summary(fit_crl)$r.squared)
  r2_crl_boxcox <- append(r2_crl_boxcox, summary(fit_crl_box)$r.squared)

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


## save r2 plots and csv
r2_df <- data.frame(labels, r2_whole, r2_whole_box_cox, r2_whole_log, r2_whole_log_cox, r2_crl, r2_crl_boxcox)
colnames(r2_df) <- c('label', 'whole_embryo', 'whole_embryo_bc', 'whole_embryo_log', 'whole_embryo_log_box', 'crl', 'crl_bc')

write.csv(r2_df, file.path(out_dir, 'r2.csv'))

boxplot_path = file.path(out_dir, 'r2.png')
boxplot(r2_df)
png(file=boxplot_path)
dev.off()




