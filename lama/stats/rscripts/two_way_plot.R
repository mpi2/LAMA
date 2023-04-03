# Title     : TODO
# Objective : TODO
# Created by: u5823099
# Created on: 21/03/2022


if (!require(grid)) install.packages(c('dplyr','factoextra', 'janitor', 'readr',
                                        'tidyverse', 'ggplot2','cluster','ggforce',
                                        'cowplot', 'grid','gridExtra', 'stringr'),
                                        repos='http://cran.us.r-project.org')

args <- commandArgs(trailingOnly = TRUE)

organ_file <- args[1];
staging_file <- args[2];
label_file <-args[3];
voxel_size <-args[4];


library(dplyr)
library(readr)
library(tidyverse)
library(janitor)
library(ggplot2)
library(factoextra)
library(cluster)
library(ggforce)
library(cowplot)
library(grid)
library(gridExtra)
library(stringr)

### set up some cool functions
variable_names <- list(
  "WT" = expression(bolditalic("Zic2")^bolditalic("+/+")),
  "HET" = expression(bolditalic("Zic2")^bolditalic("Ku/+"))
)

variable_labeller <- function(variable,value){return(variable_names[value])}

back_names <- list(
  "C3H" = "C3H/HeH",
  "C57BL6" = "C57BL6/N"
)

back_labeller <- function(variable,value){return(back_names[value])}




#This is the only way to call to make it work!
inter_names <- list(
  "WT.C3H" = expression(bolditalic("Zic2")^bolditalic("+/+")~bold("C3H/HeH")),
  "HET.C3H" = expression(bolditalic("Zic2")^bolditalic("Ku/+")~bold("C3H/HeH")),
  "WT.C57BL6" = expression(bolditalic("Zic2")^bolditalic("+/+")~bold("C57BL6/N")),
  "HET.C57BL6" = expression(bolditalic("Zic2")^bolditalic("Ku/+")~bold("C57BL6/N"))
)

inter_labeller <- function(variable,value){return(inter_names[value])}


#Functions for stat tests:

pandt_vals <- function(fit) {
  # so for some reason these are not two dimensional and only 1d

  # get estimates
  est <- fit$coefficients[fit$qr$pivot]

  # get R: see stats:::summary.lm to see how this is calculated
  p1 <- 1L:(fit$rank)
  R <- diag(chol2inv(fit$qr$qr[p1, p1, drop = FALSE]))

  # get residual sum of squares for each
  resvar <- sum(fit$residuals^2) / fit$df.residual
  # R is same for each coefficient, resvar is same within each model
  se <- sqrt(outer(R, resvar))

  tvals <- est / se
  pvals <- pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE) * 2

  return(list(pvals=pvals, tvals=tvals))
}


adj_p <- function(pvals){
  p.adjust(pvals, method='BH')
}

#converts to mm^3
scale_to_mm3 <- function(x) {
  #voxel size is 40um
  um3_conv_factor <- voxel_size^3  # To convert voxels to um3
  um3_to_mm3_conv_factor <- 1e9
  return((x * um3_conv_factor)/um3_to_mm3_conv_factor)
}

#normaliser
normalise <- function(x, na.rm = FALSE) (x/full_info$WEV)

# staging
staging_info <- read_csv(staging_file)


organ_info <- read_csv(organ_file) %>% remove_constant()

#add in factors remove empty labels
organ_info <- arrange(transform(organ_info, Genotype=factor(Genotype,levels=c("WT","HET"))))


organ_info <- arrange(transform(organ_info, Background=factor(Background,levels=c("C3H","C57BL6"))))


label_info <- read_csv(label_file)


# rename stuff
names(staging_info)[names(staging_info) == "vol"] <- "Embryo"

#get the label_info into the right format

label_info <- t(label_info)

#remove annoying first column
label_info <- label_info[,-1]

str(label_info[2,])

# merge dfs
full_info <- merge(staging_info, organ_info,
                   by = c("Embryo","Genotype", "Background"),
                   all.x = TRUE, all.y = TRUE)


full_info <- full_info %>%
  mutate_at(vars(contains('X')), scale_to_mm3) %>%
  mutate_at(vars(contains('WEV')), scale_to_mm3) %>%
  mutate_at(vars(contains('X')), normalise)

#re-order
full_info <- arrange(transform(full_info, Background=ordered(Background,levels=c("C3H","C57BL6"))))

full_info <- arrange(transform(full_info, Genotype=ordered(Genotype,levels=c("WT","HET"))))

#get proper organ_names

name_list <-c("Embryo", "Genotype", "Background", "WEV", label_info[2, ])

colnames(full_info) <- name_list

#Do full and pairwise comparisons:
g_by_e_lm <- function(y, dataset, pair_test=F){
  if (pair_test){
    library(mixlm)
    #Had to transform the data to not get residuals adding to 0?
    model <- lm(I(y * 1e6)~Genotype:Background+Genotype+Background+WEV, data=dataset)
    pair_model <- simple.glht(model,'Genotype:Background')
    #print("it's me")
    pvals <- pair_model$res$'P(>t)'
    names(pvals) <- rownames(pair_model$res)
    #gets p and t vals from model
    #results <- list(pvals=pair_model$res$'P(>t)', tvals=pair_model$res$'t value')
    detach("package:mixlm", unload=T)
    library(stats)
    return(pvals)
  }
  else{
    #dataset$Genotype <- factor(dataset$Genotype, ordered=F)
    #dataset$Background <- factor(dataset$Background, ordered=F)
    model <- lm(y~Genotype:Background+Genotype+Background+WEV, data=dataset)
    results <- pandt_vals(model)
    pvals <- results$pvals
    names(pvals) <- names(model$coefficients)
    return(pvals)
  }
}

g_by_e_pvals <- sapply(full_info[5:191], function(data)
  g_by_e_lm(data, full_info, pair_test = F)
)

g_by_e_qvals <- t(apply(g_by_e_pvals[c(2,3,5), ], 1, function(x)  p.adjust(x, method='BH')))


pairwise_pvals <- sapply(full_info[5:191], function(data)
  g_by_e_lm(data, full_info, pair_test = T)
)

pairwise_qvals <- t(apply(pairwise_pvals[c(1,2,5,6), ], 1, function(x)  p.adjust(x, method='BH')))


# function to annotate_pvals
anno_qvals <- function(qvals) {
  paste_vect <- list(length(qvals))
  #print(names(qvals))
  for (col in seq_along(qvals)) {
    #remove the stupid signed integer crap for the g_by_e_stuff
    Names <- names(qvals[col]) %>%
      str_remove_all(".L") %>%
      str_replace("WT:C57BL6", as.character("wildtype C57BL6/N")) %>%
      str_replace("WT:C3H", as.character("wildtype C3H/HeH")) %>%
      str_replace("HET:C57BL6", as.character("mutant C57BL6/N")) %>%
      str_replace("HET:C3H", as.character("mutant C3H/HeH")) %>%
      str_replace(":", " X ") %>%
      str_replace("-", " v.s. ")
    paste_vect[col] <- paste0(Names, " = ", round(qvals[col], 6))
  }
  return(paste_vect)
}

g_by_e_anno <- lapply(seq_along(g_by_e_qvals[1,]), function(i)
  anno_qvals(g_by_e_qvals[,i])
)

pair_anno <- lapply(seq_along(pairwise_qvals[1,]), function(i)
  anno_qvals(pairwise_qvals[,i])
)

###summary_grid###

#this may be stupid

rect_chooser <- function(Genotype, Background){
  ### this stupid function just makes a rectangle dataframe (i.e. left, right, top, bottom) beased on the condition per row
  xleft = vector()
  xright = vector()
  ybottom = vector()
  ytop = vector()
  for (i in (1:length(Genotype))){
    if (interaction(Genotype[i],Background[i])=="WT.C3H"){
      xleft <- append(xleft, 0)
      xright <- append(xright, 0.9)
      ybottom <- append(ybottom, 1.1)
      ytop <- append(ytop, 2)
    }
    else if (interaction(Genotype[i],Background[i])=="HET.C3H"){
      xleft <-append(xleft, 1.1)
      xright <- append(xright, 2)
      ybottom <-append(ybottom, 1.1)
      ytop <- append(ytop, 2)
    }
    else if (interaction(Genotype[i],Background[i])=="WT.C57BL6"){
      xleft <- append(xleft, 0)
      xright <- append(xright, 0.9)
      ybottom <- append(ybottom, 0)
      ytop <- append(ytop, 0.9)
    }
    else{
      xleft <- append(xleft, 1.1)
      xright <- append(xright, 2)
      ybottom <- append(ybottom, 0)
      ytop <-append(ytop, 0.9)
    }

  }
  return(data.frame(l = xleft, r = xright, b = ybottom, t = ytop))
}



# So get the rectangle values
dumb_sum <- full_info[1:3] %>%
  mutate(Rect=rect_chooser(Genotype, Background))


# now design p-val
square_plots <- lapply(seq_along(full_info[5:191]), function(i)
  ggplot()+
    geom_rect(dumb_sum, mapping=aes(xmin=Rect$l, xmax=Rect$r, ymin=Rect$b, ymax=Rect$t,
                                    fill=interaction(Genotype,Background)))+
    geom_text(dumb_sum, mapping=aes(x=(Rect$l+Rect$r)/2, y=(Rect$b+Rect$t)/2,  # add labels in center of rectangle- needs brackets to work
                                    label = sapply(interaction(Genotype,Background),
                                                   function(cond) gsub("expression","", inter_labeller(1, cond)) )), # calls inter_labeller on just the one row and
              # removes the expression
              parse=T)+
    {if(pairwise_qvals[1, i] < 0.05 ) geom_segment(aes(x = -0.1, xend = 2.1, y = 2.1, yend = 2.1), # wt veh vs het veh
                                                   arrow = arrow(type = "closed",ends = "both"))}+
    {if(pairwise_qvals[2, i] < 0.05 ) geom_segment(aes(x = -0.1, xend = -0.1, y = -0.1, yend = 2.1), # wt veh vs wt eth
                                                   arrow = arrow(type = "closed",ends = "both"))}+
    {if(pairwise_qvals[3, i] < 0.05 ) geom_segment(aes(x = 2.1, xend = 2.1, y = -0.1, yend = 2.1), # het veh vs het eth
                                                   arrow = arrow(type = "closed",ends = "both"))}+
    {if(pairwise_qvals[4, i] < 0.05) geom_segment(aes(x = -0.1, xend = 2.1, y = -0.1, yend = -0.1), # wt eth vs het eth
                                                  arrow = arrow(type = "closed",ends = "both"))}+
    {if(g_by_e_qvals[1, i] < 0.05) geom_segment(aes(x = 1, xend = 1, y = 0, yend = 2, color = "red"))}+
    {if(g_by_e_qvals[1, i] < 0.05) geom_curve(aes(x = 0.45, xend = 1.55, y = 1, yend = 1),
                                              curvature = 0.25, arrow = arrow(type = "closed", ends = "both"))}+

    {if(g_by_e_qvals[2, i] < 0.05) geom_segment(aes(x = 0, xend = 2, y = 1, yend = 1, color = "red"))}+
    {if(g_by_e_qvals[2, i] < 0.05) geom_curve(aes(x = 1, xend = 1, y = 0.45, yend = 1.55),
                                              curvature = 0.25, arrow = arrow(type = "closed",ends = "both"))}+

    {if(g_by_e_qvals[3, i] < 0.05) geom_ellipse(aes(x0=1.55,y0=0.45,a=0.5,b=0.5, angle=0, colour="red"))}+
    #ggtitle(stringr::str_to_title(gsub("_", " ", names(full_info[5:191])[[i]])))+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+#,plot.title = element_text(hjust = 0.5, size = 40),)+
    scale_fill_manual(name="", labels=unlist(inter_names), values=c("#F8766D", "#00BFC4","#FFA500","#C77CFF"))+
    labs(tag = "A")

)# scales to colours we want

### Draw g_by_e_boxplots for results
box_plots <- lapply(seq_along(full_info[5:191]), function(data)
  ggplot(full_info,aes(x=WEV,y=full_info[5:191][[data]],colour=interaction(Genotype, Background),label=Embryo))+
    geom_boxplot()+
    geom_point()+
    #geom_smooth(method = "lm")+
    facet_grid(~interaction(Genotype,Background), labeller = inter_labeller, scales = "free_x")+
    geom_text(aes(label=ifelse(data > quantile(data, 0.975),as.character(Embryo),'' )), hjust=0, vjust=0)+
    geom_text(aes(label=ifelse(data < quantile(data, 0.025),as.character(Embryo),'' )), hjust=0, vjust=0)+
    scale_color_manual(name="", labels=unlist(inter_names), values=c("#F8766D", "#00BFC4","#FFA500","#C77CFF"))+
    ylab(expression("Organ Volume Normalised To WEV"))+
    xlab(expression("Whole Embryo Volume (WEV)"))+
    ggtitle(paste(g_by_e_anno[data][[1]][[1]], " ",
                  g_by_e_anno[data][[1]][[2]]) , subtitle = g_by_e_anno[data][[1]][[3]])+
    theme(plot.title = element_text(hjust = 0.5, size = 10), plot.subtitle=element_text(size=10, hjust=0.5),
          strip.text.x = element_text(hjust = 0.5, size = 6), legend.position = "none" )+
    labs(tag = "B"))


## genotype pairwise comparisons
geno_line_plots <- lapply(seq_along(full_info[5:191]), function(data)
  ggplot(full_info,aes(x=WEV,y=full_info[5:191][[data]],colour=interaction(Genotype, Background),label=Embryo))+
    #geom_boxplot()+
    geom_point()+
    geom_smooth(method = "lm")+
    facet_grid(~Background, labeller = back_labeller, scales = "free_x")+
    geom_text(aes(label=ifelse(data > quantile(data, 0.975),as.character(Embryo),'' )), hjust=0, vjust=0)+
    geom_text(aes(label=ifelse(data < quantile(data, 0.025),as.character(Embryo),'' )), hjust=0, vjust=0)+
    scale_color_manual(name="", labels=unlist(inter_names), values=c("#F8766D", "#00BFC4","#FFA500","#C77CFF"))+
    ylab(expression("Organ Volume Normalised To WEV"))+
    xlab(expression("Whole Embryo Volume (WEV)"))+
    labs(tag = "C")+
    ggtitle(pair_anno[data][[1]][[1]], subtitle = pair_anno[data][[1]][[4]])+
    theme(plot.title = element_text(hjust = 0.5, size = 10), plot.subtitle=element_text(size=10, hjust=0.5),legend.position = "none"))

# annotate facets with p-vals
# geno_tagged <-lapply(seq_along(geno_line_plots), function(p)
#   tag_facet(geno_line_plots[p][[1]],
#             tag_pool = c(
#               as.character(pair_anno[p][[1]][[1]]),
#               as.character(pair_anno[p][[1]][[4]])),
#             open="", close="",
#             fontface = 4,
#             size=4)
# )

### Background pairwise p-vals
treat_line_plots <- lapply(seq_along(full_info[5:191]), function(data)
  ggplot(full_info,aes(x=WEV,y=full_info[5:191][[data]],colour=interaction(Genotype, Background),label=Embryo))+
    #geom_boxplot()+
    geom_point()+
    geom_smooth(method = "lm")+
    facet_grid(~Genotype, labeller = variable_labeller, scales = "free_x")+
    geom_text(aes(label=ifelse(data > quantile(data, 0.975),as.character(Embryo),'' )), hjust=0, vjust=0)+
    geom_text(aes(label=ifelse(data < quantile(data, 0.025),as.character(Embryo),'' )), hjust=0, vjust=0)+
    scale_color_manual(name="", labels=unlist(inter_names), values=c("#F8766D", "#00BFC4","#FFA500","#C77CFF"))+
    ylab(expression("Organ Volume Normalised To WEV"))+
    xlab(expression("Whole Embryo Volume (WEV)"))+
    labs(tag = "D")+
    ggtitle(pair_anno[data][[1]][[2]], subtitle = pair_anno[data][[1]][[3]])+
    theme(plot.title = element_text(hjust = 0.5, size = 10), plot.subtitle=element_text(size=10, hjust=0.5),legend.position = "none"))


#### PCA analysis ####

# basic PCA
PCAs <- lapply(full_info[5:191], function(data)
  prcomp(data.frame(org=data, WEV=full_info[4]), scale=T)
)


groups <- interaction(full_info$Genotype, full_info$Background)



PCA_plots <- lapply(PCAs,function(data)
  fviz_pca_ind(data,
               geom = "point",
               legend.title = "Groups",
               pallete =c("#F8766D", "#00BFC4","#FFA500","#C77CFF"),
               habillage = interaction(full_info$Genotype, full_info$Background),
               addEllipses = T,
               repel = TRUE,
               title="PCA analysis")+
    scale_fill_manual(name="", labels=unlist(inter_names), values=c("#F8766D", "#00BFC4","#FFA500","#C77CFF"))+
    scale_color_manual(name="", labels=unlist(inter_names), values=c("#F8766D", "#00BFC4","#FFA500","#C77CFF"))+
    theme(legend.position = "none")+
    labs(tag = "E")
)

#### Play around with clustering ####

df <- (data.frame(org=full_info[7], WEV=full_info[4]))

# Set up test data-frame to develop clustering. Scale as organ volumes are fractions while whole embryo volumes are massive


my_list <- list("1"="1", "2"="2", "3" = "3", "4"="4")

dfs <- lapply(full_info[5:191], function(data)
  scale(data.frame(org=data, WEV=full_info[4])))

km_plots <- lapply(dfs, function(data){
  fviz_cluster(kmeans(data, 4, nstart = 25), data = data, geom = "point", repel = T, main="Kmeans Cluster")})


for(i in seq_along(km_plots)){
  km_plots[[i]]$Group <- interaction(full_info$Genotype, full_info$Background)
  str(km_plots[[i]]$Group)
}


km_plots_lab <- lapply(km_plots, function(plot)
  plot + geom_point(aes(colour = plot$Group))+
    scale_fill_manual(values = c("white", "white","white","white"))+
    scale_color_manual(labels=c(my_list,
                                inter_names),
                       values = c("black", "black","black","black", "#C77CFF", "#00BFC4","#FFA500", "#F8766D"))+
    theme(legend.position = "none")+
    labs(tag = "F") )

#Dumb fix but it will do
for(i in seq_along(km_plots)){
  km_plots[[i]]$data$Group <- interaction(full_info$Genotype, full_info$Background)
  str(km_plots[[i]]$Group)
}

cluster.summaries <- lapply(km_plots, function(pl)
  pl$data %>%
    group_by(Group, cluster) %>%
    summarise(n=n()) %>%
    group_by(Group) %>%
    mutate(nGroup = sum(n), Percentage = n / sum(n) * 100)
)



#Make pdf:

pdf("g_by_back_plots_v1.pdf", onefile = T, paper="a4r", width=13, height=10)
for(i in seq_along(square_plots)){

  #plot rect summary in square_1
  p1 <- as_grob(square_plots[[i]])
  # plot boxplots in square_2 if g_by_e signif
  if (any(g_by_e_qvals[,i] < 0.05)){
    p2 <- as_grob(box_plots[[i]])
  }
  else {p2 <-as_grob(geom_blank())}
  # plot geno line plots in square 3 if g-pairwise are significant
  if (any(pairwise_qvals[c(1,4), i] < 0.05)){
    p3 <- as_grob(geno_line_plots[[i]])
  }
  else{p3<-as_grob(geom_blank())}
  # same thing for treat line plots in square4
  if (any(pairwise_qvals[c(2,3), i] < 0.05)){
    p4 <- as_grob(treat_line_plots[[i]])
  }
  else{ p4 <-as_grob(geom_blank())}

  if (any(cluster.summaries[[i]]$Percentage > 60.5) ) {
    p5 <- as_grob(PCA_plots[[i]])

    p6 <- as_grob(km_plots_lab[[i]])
  }
  else{
    p5 <-as_grob(geom_blank())
    p6 <-as_grob(geom_blank())
  }

  grid.arrange(grobs=list(p1, p2, p3, p4, p5, p6), ncol=3, nrow=2,
               top = textGrob(stringr::str_to_title(gsub("_", " ", names(full_info[5:191])[[i]])),  gp=gpar(fontsize=28,font=8)))
}

dev.off()








