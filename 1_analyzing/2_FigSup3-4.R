rm(list=ls())
InputFolder <- ".../0_file/"
#import functions used for statistics
source(paste0(InputFolder,"0_functions.R"))
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, scales, gtools, lme4, lmerTest, ggsci, rcarbon, ppcor,
               psych, agricolae, cowplot, ggpubr)
df.main <- read.csv(paste0(InputFolder, "main.csv"))
# 1 - data cleaning and filtering ###############################################
df.clean <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                               daytime_VPD < 2,!(PET %in% c(99,12,7,9))) %>% drop_na(Amax2000) 

df.clean <-transform(df.clean, bin_fAPAR = cut(df.clean$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean <-transform(df.clean, bin_T = cut(df.clean$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))

df.clean.cor <- data.frame(FLUXNET = df.clean$Amax2000,
                           LagPAR = df.clean$LagPAR,
                           LagTair = df.clean$LagT,
                           VPD = df.clean$daytime_VPD,
                           LagVPD = df.clean$LagVPD,
                           Alpha = df.clean$daytime_Alpha,
                           Tair = df.clean$daytime_Tair,
                           fAPAR = df.clean$fAPAR,
                           bin_fAPAR = df.clean$bin_fAPAR, 
                           bin_T = df.clean$bin_T,
                           Site = as.factor(df.clean$Site),
                           PFT = as.factor(df.clean$PET))
# 2 - extract info from each bin ###############################################
BinInfo <- data.frame(ddply(df.clean.cor, .(bin_fAPAR,  bin_T), bin.info))

limit_plot <- readRDS(paste0(InputFolder, "fPAR_range.rds"))
colrev <- rev(heat.colors(999))

FigSup3a <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = LagTmean))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,30), oob=squish,
                       name=expression('(\u00B0C)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "a") +
  ggtitle(bquote(bar(T[air])~'(mean)')) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup3b <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = LagTsd))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,6), oob=squish,
                       name=expression('(\u00B0C)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "b") +
  ggtitle(bquote(bar(T[air])~'(SD)')) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup3c <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = LagTrange))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,12.5), oob=squish,
                       name=expression('(\u00B0C)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "c") +
  ggtitle(bquote(bar(T[air])~'(range)')) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup3d <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = LagPARmean))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(145,450), oob=squish,
                       name=bquote("("*mu*"mol"~m^-2~s^-1*")")) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "d") +
  ggtitle(bquote(bar(PPFD)~'(mean)')) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup3e <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = LagVPDmean))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,2), oob=squish,
                       name= "(kPa)") +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "e") +
  ggtitle(bquote(bar(VPD)~'(mean)')) +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup3f <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = VPDmean))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,2), oob=squish,
                       name= "(kPa)") +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "f") +
  ggtitle("VPD (mean)") +
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup3 <- plot_grid(FigSup3a, FigSup3b, FigSup3c, 
                     FigSup3d, FigSup3e, FigSup3f, 
                     align = 'hv', ncol = 3, nrow = 2)
# 3 - extract PFT info from each bin ###############################################
FigSup4a <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = DBF*100))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,100), oob=squish,
                       name=expression('(%)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "a") +
  ggtitle(expression('DBF'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup4b <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = EBF*100))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,100), oob=squish,
                       name=expression('(%)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "b") +
  ggtitle(expression('EBF'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup4c <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = ENF*100))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,100), oob=squish,
                       name=expression('(%)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "c") +
  ggtitle(expression('ENF'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup4d <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = GRA*100))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,100), oob=squish,
                       name=expression('(%)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "d") +
  ggtitle(expression('GRA'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup4e <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = MF*100))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,100), oob=squish,
                       name=expression('(%)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "e") +
  ggtitle(expression('MF'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))

FigSup4f <- ggplot(data = BinInfo, aes(bin_fAPAR, bin_T, fill = WET*100))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradientn(colours = colrev, limit = c(0,100), oob=squish,
                       name=expression('(%)')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  labs(tag = "f") +
  ggtitle(expression('WET'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(
    barwidth = 1, barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"))
