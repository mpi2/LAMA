
df <- read.csv('/home/neil/bit/LAMA_results/E14.5/paper_runs/mutant_runs/280618_analysed_lines/nras/output/stats/organvolumes/specimen_calls/20170125_NRAS_E14.5_4.2b_HOM_XX_REC_scaled_4.7297_pixel_13.9999.nrrd_inverted_organ_volumes_LM_FDR5%.csv')

out <- p.adjust(df$p, method='BH')