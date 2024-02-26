rm(list=ls())
InputFolder <- ".../0_file/"
#import functions used for statistics
source(paste0(InputFolder,"0_functions.R"))
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, scales, gtools, lme4, lmerTest, ggsci, rcarbon, ppcor,
               psych, agricolae, cowplot, gridExtra, ggpubr)

df.main <- read.csv(paste0(InputFolder, "main.csv"))
# 1 - Analysis for each PFT ----------------------------------------------
df.clean <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                               daytime_VPD < 2,!(PET %in% c(99,7,9))) %>% drop_na(Amax2000)

df.clean <-transform(df.clean, bin_fAPAR = cut(df.clean$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean <-transform(df.clean, bin_T = cut(df.clean$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))

df.clean.cor <- data.frame(bin_fAPAR = df.clean$bin_fAPAR,
                           bin_T = df.clean$bin_T,
                           fAPAR = df.clean$fAPAR,
                           Tair = df.clean$daytime_Tair,
                           Site = as.factor(df.clean$Site))

list.frames <- replicate(7, data.frame())
df.pft <- data.frame()

for (i in 1: length(unique(df.clean$PET))) {
  
  PFT <- unique(df.clean$PET)[i]
  df.clean.pft <- df.clean[(df.clean$PET == PFT),] %>% filter(fAPAR < 0.98)
  
  if (PFT == 12) {Group = "CRO"}
  if (PFT == 4) {Group = "DBF"}
  if (PFT == 2) {Group = "EBF"}
  if (PFT == 1) {Group = "ENF"}
  if (PFT == 5) {Group = "MF"} 
  if (PFT == 10) {Group = "GRA"} 
  if (PFT == 11) {Group = "WET"} 
  
  df.clean.pft <-transform(df.clean.pft, bin_fAPAR = cut(df.clean.pft$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
  df.clean.pft <-transform(df.clean.pft, bin_T = cut(df.clean.pft$daytime_Tair, c(seq(0, 38, by = 1), Inf), include.lowest = TRUE))
  
  df.clean.pft.cor <- data.frame(FLUXNET = df.clean.pft$Amax2000, 
                                         LagPAR = df.clean.pft$LagPAR,
                                         LagTair = df.clean.pft$LagT,
                                         bin_fAPAR = df.clean.pft$bin_fAPAR, 
                                         bin_T = df.clean.pft$bin_T,
                                         PFT = as.factor(df.clean.pft$PET),
                                         Site = as.factor(df.clean.pft$Site))
  
  LMM.run.pft <- data.frame(ddply(df.clean.pft.cor, .(bin_fAPAR,bin_T), LMM_PFT))
  LMM.run.pft$Group = Group
  LMM.run.pft$sig.flag <- cut(LMM.run.pft$pValue,breaks=c(-Inf, 0.05, Inf), label=c("\u2022", "")) 
  
  list.frames[[i]] <- LMM.run.pft
  df.pft <- rbind(LMM.run.pft, df.pft)
}

limit_fAPAR <- sort(unique(df.clean.cor$bin_fAPAR))
limit_T <- sort(unique(df.clean.cor$bin_T))

df.clean.base <- data.frame(bin_fAPAR = df.clean.cor$bin_fAPAR, 
                            bin_T = df.clean.cor$bin_T,
                            fAPAR = df.clean.cor$fAPAR, 
                            Tair = df.clean.cor$Tair, 
                            reg.coef = 0, 
                            pValue = 1)
t.mf <- list.frames[[1]]
t.enf <- list.frames[[2]]
t.gra <- list.frames[[3]]
t.ebf <- list.frames[[4]]
t.wet <- list.frames[[5]]
t.dbf <- list.frames[[6]]
t.cro <- list.frames[[7]]

# CRO
cro.d <- round(nrow(t.cro[t.cro$reg.coef>0,])/nrow(t.cro),2)*100
fig.ext.a <- ggplot(data = t.cro, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T,
                    breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]","(35,36]"),
                   labels = c("0","5","10","15","20","25","30","35")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "a") +
  ggtitle(paste0("CRO (", cro.d, "%)"))+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_colourbar(title.position = "right",
                                barwidth = 1, barheight = 10,
                                title.hjust = 0.5,frame.colour = "black",
                                ticks.colour = "black"))
# DBF
dbf.d <- round(nrow(t.dbf[t.dbf$reg.coef>0,])/nrow(t.dbf)*100,0)
fig.ext.b <- ggplot(data = t.dbf, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T[1:31],
                   breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "b") +
  ggtitle(paste0("DBF (", dbf.d, "%)"))+
  theme_bw()+ 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# EBF
ebf.d <- round(nrow(t.ebf[t.ebf$reg.coef>0,])/nrow(t.ebf)*100,0)
fig.ext.c <- ggplot(data = t.ebf, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T[1:31],
                   breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]","(35,36]"),
                   labels = c("0","5","10","15","20","25","30","35")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "c") +
  ggtitle(paste0("EBF (", ebf.d, "%)"))+
  theme_bw()+ 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# ENF
enf.d <- round(nrow(t.enf[t.enf$reg.coef>0,])/nrow(t.enf)*100,0)
fig.ext.d <- ggplot(data = t.enf, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR[1:31],
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T,
                   breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]","(35,36]"),
                   labels = c("0","5","10","15","20","25","30","35")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "d") +
  ggtitle(paste0("ENF (", enf.d, "%)"))+
  theme_bw()+ 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# GRA
gra.d <- round(nrow(t.gra[t.gra$reg.coef>0,])/nrow(t.gra)*100,0)
fig.ext.e <- ggplot(data = t.gra, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR[1:31],
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T[1:31],
                   breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]","(35,36]"),
                   labels = c("0","5","10","15","20","25","30","35")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "e") +
  ggtitle(paste0("GRA (", gra.d, "%)"))+
  theme_bw()+ 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# MF
mf.d <- round(nrow(t.mf[t.mf$reg.coef>0,])/nrow(t.mf)*100,0)
fig.ext.f <- ggplot(data = t.mf, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T[1:31],
                   breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]","(35,36]"),
                   labels = c("0","5","10","15","20","25","30","35")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "f") +
  ggtitle(paste0("MF (", mf.d, "%)"))+
  theme_bw()+ 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# WET
wet.d <- round(nrow(t.wet[t.wet$reg.coef>0,])/nrow(t.wet)*100,0)
fig.ext.g <- ggplot(data = t.wet, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", size = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-5,5), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_fAPAR,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]",
                              "(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(limits = limit_T[1:31],
                   breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]",
                              "(25,26]","(30,31]","(35,36]"),
                   labels = c("0","5","10","15","20","25","30","35")) +
  geom_text(aes(label=sig.flag), color="black", size= 4) + 
  labs(tag = "g") +
  ggtitle(paste0("WET (", wet.d, "%)"))+
  theme_bw()+ 
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(angle = -90),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

