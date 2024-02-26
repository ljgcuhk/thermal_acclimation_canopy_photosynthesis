rm(list=ls())
InputFolder <- ".../0_file/"
#import functions used for statistics
source(paste0(InputFolder,"0_functions.R"))
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, lme4, lmerTest, scales, ggpubr)
df.main <- read.csv(paste0(InputFolder, "main.csv"))
# 0 - data clean ----------------------------------------------
df.clean <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                               daytime_VPD < 2,!(PET %in% c(99,12,7,9))) %>% drop_na(Amax2000) # 149,774 samples are left

# group the data set by fAPAR (range 0.3-1 by 0.02) and by Tair (range 0-35 by 1)
df.clean <-transform(df.clean, bin_fAPAR = cut(df.clean$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean <-transform(df.clean, bin_T = cut(df.clean$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))

# create a cleaner dataset for further correlation analysis
df.clean.cor <- data.frame(FLUXNET = df.clean$Amax2000, 
                           LagTair = df.clean$LagT,
                           bin_fAPAR = df.clean$bin_fAPAR, 
                           bin_T = df.clean$bin_T,
                           Site = as.factor(df.clean$Site))

# 1 - Analysis---------------------------------------------------
LMM.run <- data.frame(ddply(df.clean.cor, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run$sig.flag <- cut(LMM.run$pValue, breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

limit_plot <- readRDS(paste0(InputFolder, "fPAR_range.rds"))
colrev <- rev(heat.colors(999))
LMM.run$marginal.r <- sqrt(LMM.run$marginal.r2)
LMM.run$conditional.r <- sqrt(LMM.run$conditional.r)

FigSup5a <- ggplot(data = LMM.run, aes(bin_fAPAR, bin_T, fill = N.sample))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,500), oob=squish) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                  breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "a") +
  ggtitle("Sample number")+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup5b <- ggplot(data = LMM.run, aes(bin_fAPAR, bin_T, fill = marginal.r))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,0.61), oob=squish) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
    breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "b") +
  ggtitle(expression(paste("Marginal ",italic("r"))))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup5c <- ggplot(data = LMM.run, aes(bin_fAPAR, bin_T, fill = conditional.r))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0.3,1), oob=squish) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
    breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "c") +
  ggtitle(expression(paste("Conditional ",italic("r"))))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

