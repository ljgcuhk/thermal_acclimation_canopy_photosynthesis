rm(list=ls())
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, scales, gtools, lme4, lmerTest, ggsci, rcarbon, ppcor,
               psych, agricolae, cowplot)
# 1 - read data----------------------------------------------------------------
InputFolder <- ".../0_file/"
#import functions used for statistics
source(paste0(InputFolder,"0_functions.R"))
df.main <- read.csv(paste0(InputFolder, "main.csv"))

# 1 - Analysis: Amax (not standardized to PPFD = 2000)---------------------------------------------
df.clean.Amax <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                                    daytime_VPD < 2,!(PET %in% c(99,12,7,9)))%>% drop_na(Amax) 

df.clean.Amax <-transform(df.clean.Amax, bin_fAPAR = cut(df.clean.Amax$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean.Amax <-transform(df.clean.Amax, bin_T = cut(df.clean.Amax$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))

df.clean.Amax.cor <- data.frame(FLUXNET = df.clean.Amax$Amax, 
                           LagTair = df.clean.Amax$LagT,
                           bin_fAPAR = df.clean.Amax$bin_fAPAR, 
                           bin_T = df.clean.Amax$bin_T,
                           PFT = as.factor(df.clean.Amax$PET),
                           Site = as.factor(df.clean.Amax$Site))

LMM.run.Amax <- data.frame(ddply(df.clean.Amax.cor, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run.Amax$sig.flag <- cut(LMM.run.Amax$pValue, breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

# positive (n > 20)  = 0.8732694
print(nrow(LMM.run.Amax[LMM.run.Amax$reg.coef>0,])/(nrow(LMM.run.Amax[LMM.run.Amax$reg.coef>0,])+nrow(LMM.run.Amax[LMM.run.Amax$reg.coef<0,])))

limit_plot <- readRDS(paste0(InputFolder, "fPAR_range.rds"))

FigExt2c <- ggplot(data = LMM.run.Amax, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-4,4), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (°C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                  breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  ggtitle(expression('A'[max]~'~'~bar(T[air])*' + (1|Site) (87%)'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(label.position = "right",
                                title.position = "right",
                                barwidth = 0.8, barheight = 10,
                                title.hjust = 0.5,frame.colour = "black",
                                ticks.colour = "black"))

# 2 - Analysis: Amax2000 with PFT for random effect-----------------------------
df.clean.Amax2000 <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                                    daytime_VPD < 2,!(PET %in% c(99,12,7,9)))%>% drop_na(Amax2000) 

df.clean.Amax2000 <-transform(df.clean.Amax2000, bin_fAPAR = cut(df.clean.Amax2000$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean.Amax2000 <-transform(df.clean.Amax2000, bin_T = cut(df.clean.Amax2000$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))

df.clean.Amax2000.cor <- data.frame(FLUXNET = df.clean.Amax2000$Amax2000, 
                                LagTair = df.clean.Amax2000$LagT,
                                bin_fAPAR = df.clean.Amax2000$bin_fAPAR, 
                                bin_T = df.clean.Amax2000$bin_T,
                                PFT = as.factor(df.clean.Amax2000$PET),
                                Site = as.factor(df.clean.Amax2000$Site))

LMM.run.Amax2000 <- data.frame(ddply(df.clean.Amax2000.cor, .(bin_fAPAR,  bin_T), LMM_ALL_PFT))
LMM.run.Amax2000$sig.flag <- cut(LMM.run.Amax2000$pValue, breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

# positive (n > 20)  = 0.85836
print(nrow(LMM.run.Amax2000[LMM.run.Amax2000$reg.coef>0,])/(nrow(LMM.run.Amax2000[LMM.run.Amax2000$reg.coef>0,])+nrow(LMM.run.Amax2000[LMM.run.Amax2000$reg.coef<0,])))

FigExt2d <- ggplot(data = LMM.run.Amax2000, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-2,2), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (°C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  ggtitle(expression('A'[max*",2000"]~'~'~bar(T[air])*' + (1|PFT) (86%)'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(label.position = "right",
                                title.position = "right",
                                barwidth = 0.8, barheight = 10,
                                title.hjust = 0.5,frame.colour = "black",
                                ticks.colour = "black"))


# 3 - Analysis: Amax2000 ~ GrowthT + GrowthPPFD + (1|Site) ----------------------
df.clean.pcoef.cor <- data.frame(FLUXNET = df.clean.Amax2000$Amax2000, 
                                LagTair = df.clean.Amax2000$LagT,
                                LagPAR = df.clean.Amax2000$LagPAR,
                                bin_fAPAR = df.clean.Amax2000$bin_fAPAR, 
                                bin_T = df.clean.Amax2000$bin_T,
                                Site = as.factor(df.clean.Amax2000$Site))
# derive partial coef for growth Tair
LMM.run.pcoef <- data.frame(ddply(df.clean.pcoef.cor, .(bin_fAPAR,  bin_T), LMM_T_PPFD))
LMM.run.pcoef$sig.flag <- cut(LMM.run.pcoef$pValue.temp, breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

# positive (n > 20)  = 0.8892439
print(nrow(LMM.run.pcoef[LMM.run.pcoef$reg.coef.temp>0,])/(nrow(LMM.run.pcoef[LMM.run.pcoef$reg.coef.temp>0,])+nrow(LMM.run.pcoef[LMM.run.pcoef$reg.coef.temp<0,])))

FigExt2a <- ggplot(data = LMM.run.pcoef, aes(bin_fAPAR, bin_T, fill = reg.coef.temp))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-2,2), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (°C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  ggtitle(expression('A'[max*",2000"]~'~'~bar(T[air])~'+'~bar(PPFD)~'+ (1|Site) (89%)'))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(label.position = "right",
                                title.position = "right",
                                barwidth = 0.8, barheight = 10,
                                title.hjust = 0.5,frame.colour = "black",
                                ticks.colour = "black"))


# 4 - Analysis: Partial correlation (controling for PPFD)---------------------------------------------
df.clean.pcor.cor <- data.frame(FLUXNET = df.clean.Amax2000$Amax2000, 
                                   LagPAR = df.clean.Amax2000$LagPAR,
                                   LagTair = df.clean.Amax2000$LagT,
                                   LagVPD = df.clean.Amax2000$LagVPD,
                                   bin_fAPAR = df.clean.Amax2000$bin_fAPAR, 
                                   bin_T = df.clean.Amax2000$bin_T,
                                   PFT = as.factor(df.clean.Amax2000$PET),
                                   Site = as.factor(df.clean.Amax2000$Site))

LMM.run.pcor.T <- data.frame(ddply(df.clean.pcor.cor, .(bin_fAPAR,  bin_T), Pcor.fit.T))
LMM.run.pcor.T$sig.flag <- cut(LMM.run.pcor.T$pValue,breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

# positive (n > 20)  = 0.85623
print(nrow(LMM.run.pcor.T[LMM.run.pcor.T$pcorr>0,])/
        (nrow(LMM.run.pcor.T[LMM.run.pcor.T$pcorr>0,])+
           nrow(LMM.run.pcor.T[LMM.run.pcor.T$pcorr<0,])))

FigExt2b <- ggplot(data = LMM.run.pcor.T, aes(bin_fAPAR, bin_T, fill = pcorr))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-0.4000001,0.600001), space = "Lab", oob=squish,
                       name=expression(paste("Partial ",italic("r")))) +
  ylab(expression(paste(T["air"], " (°C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                    breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  ggtitle(expression(paste('A'[max*",2000"]~'~'~bar(T[air]),"|",bar(PPFD), " (86%)")))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="right",
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(label.position = "right",
                                title.position = "right",
                                barwidth = 0.8, barheight = 10,
                                title.hjust = 0.5,frame.colour = "black",
                                ticks.colour = "black"))

