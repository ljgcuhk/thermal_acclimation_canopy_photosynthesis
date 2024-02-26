rm(list=ls())
InputFolder <- ".../0_file/"
#import functions used for statistics
source(paste0(InputFolder,"0_functions.R"))
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, scales, gtools, lme4, lmerTest, ggsci, rcarbon, ppcor,
               psych, agricolae, cowplot, ggmap, maps, RColorBrewer, ggpubr)

df.main <- read.csv(paste0(InputFolder, "main.csv"))

# 0 - data clean ----------------------------------------------
# fAPAR > 0.3 and T > 0 for filtering growing season
# Alpha = AET/PET > 0.7 for filtering drought
# PET == 7/9/12/99 for filtering drylands/cropland/snow sites
# beta is ecosystem-scale Amax. Amax2000 is beta standardized to PPFD = 2000.
# beta data has been quality controlled.
df.clean <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                               daytime_VPD < 2,!(PET %in% c(99,12,7,9))) %>% drop_na(Amax2000) 
# 149,774 samples are left
# group the data set by fAPAR (range 0.3-1 by 0.02) and by Tair (range 0-35 by 1)
df.clean <-transform(df.clean, bin_fAPAR = cut(df.clean$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean <-transform(df.clean, bin_T = cut(df.clean$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))
# create a cleaner dataset for further correlation analysis
df.clean.cor <- data.frame(FLUXNET = df.clean$Amax2000, 
                           LagTair = df.clean$LagT,
                           bin_fAPAR = df.clean$bin_fAPAR, 
                           bin_T = df.clean$bin_T,
                           PFT = as.factor(df.clean$PET),
                           Site = as.factor(df.clean$Site))

# 1 - Analysis: global evidence----------------------------------------------
# linear mixed effect model (LMM) for regression coef, pvalue, marginal r2 and conditional r2
# Amax vs LagTair with random effect by Site
# For cross-site analysis, we only report the result with N > 20
LMM.run <- data.frame(ddply(df.clean.cor, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run$sig.flag <- cut(LMM.run$pValue, breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

# check the percentage of positive (including significant) and negative coefs
print(nrow(LMM.run[LMM.run$reg.coef>0,])) #  slope > 0: 819
tmp <- LMM.run[!(LMM.run$reg.coef<0),] 
print(nrow(tmp[tmp$sig.flag == "\u2022",])) #  slope > 0 & significant: 541
print(nrow(LMM.run[LMM.run$reg.coef<0,])) # slope < 0: 120
# positive (n > 20)  = 0.8722045
print(nrow(LMM.run[LMM.run$reg.coef>0,])/(nrow(LMM.run[LMM.run$reg.coef>0,])+nrow(LMM.run[LMM.run$reg.coef<0,])))
# significant (n > 20)  = 0.6605617
print(nrow(tmp[tmp$sig.flag == "\u2022",])/nrow(LMM.run[LMM.run$reg.coef>0,]))
# acclimation rate for the cells showing positive values
mean(LMM.run[LMM.run$reg.coef>0,]$reg.coef) #0.5687891
sd(LMM.run[LMM.run$reg.coef>0,]$reg.coef) #0.3269472
# acclimation rate for all cells
mean(LMM.run$reg.coef) #0.4075306
sd(LMM.run$reg.coef) #0.6085669
# limit_plot <- unique(LMM.run$bin_fAPAR)
limit_plot <- readRDS(paste0(InputFolder, "fPAR_range.rds"))

fig2a <- ggplot(data = LMM.run, aes(bin_fAPAR, bin_T, fill = reg.coef))+
  geom_tile(color = "gray45", linewidth = 0.5)+
  scale_fill_gradient2(low = "#3C5488B2", high = "#DC0000B2", mid = "white", 
                       midpoint = 0, limit = c(-2,2), space = "Lab", oob=squish,
                       name=expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')')) +
  ylab(expression(paste(T["air"], " (\u00B0C)")))+
  xlab("fAPAR") + 
  scale_x_discrete(limits = limit_plot,
                   breaks = c("[0.3,0.32]","(0.4,0.42]","(0.5,0.52]","(0.6,0.62]","(0.7,0.72]","(0.8,0.82]","(0.9,0.92]","(0.98,1]"),
                   labels = c("0.3","0.4","0.5","0.6", "0.7","0.8","0.9","1")) +
  scale_y_discrete(breaks = c("[0,1]","(5,6]","(10,11]","(15,16]","(20,21]","(25,26]","(30,31]"),
                   labels = c("0","5","10","15","20","25","30")) +
  geom_text(aes(label=sig.flag), color="black", size= 5) + 
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="left",
        legend.title = element_text(angle = -270),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_colourbar(label.position = "left",
                                title.position = "left",
                                barwidth = 0.8, barheight = 10,
                                title.hjust = 0.5,frame.colour = "black",
                                ticks.colour = "black"))

# 2 - Analysis: PFT-sepecific gammaT----------------------------------------------
# include CRO for PFT analysis
df.clean.thermal <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT >0,
                                       LagAlpha >= 0.7, daytime_VPD < 2, 
                                       !(PET %in% c(99,7,9))) %>% drop_na(Amax2000) 

df.clean.thermal <-transform(df.clean.thermal, bin_fAPAR = cut(df.clean.thermal$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean.thermal <-transform(df.clean.thermal, bin_T = cut(df.clean.thermal$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))

list.frames.thermal <- replicate(7, data.frame())
df.pft.thermal <- data.frame()
# run for each PFT
for (i in 1: length(unique(df.clean.thermal$PET))) {
  PFT <- unique(df.clean.thermal$PET)[i]
  df.clean.thermal.pft <- df.clean.thermal[(df.clean.thermal$PET == PFT),]
  
  if (PFT == 12) {Group = "CRO"}
  if (PFT == 4) {Group = "DBF"}
  if (PFT == 2) {Group = "EBF"}
  if (PFT == 1) {Group = "ENF"}
  if (PFT == 5) {Group = "MF"} 
  if (PFT == 10) {Group = "GRA"} 
  if (PFT == 11) {Group = "WET"} 
  
  df.clean.thermal.pft <-transform(df.clean.thermal.pft, bin_fAPAR = cut(df.clean.thermal.pft$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
  df.clean.thermal.pft <-transform(df.clean.thermal.pft, bin_T = cut(df.clean.thermal.pft$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))
  
  df.clean.thermal.pft.cor <- data.frame(FLUXNET = df.clean.thermal.pft$Amax2000, 
                                         LagTair = df.clean.thermal.pft$LagT,
                                         bin_fAPAR = df.clean.thermal.pft$bin_fAPAR, 
                                         bin_T = df.clean.thermal.pft$bin_T,
                                         PFT = as.factor(df.clean.thermal.pft$PET))
  LM.run.pft <- data.frame(ddply(df.clean.thermal.pft.cor, .(bin_fAPAR,bin_T), LM_PFT))
  LM.run.pft$Group = Group
  LM.run.pft$sig.flag <- cut(LM.run.pft$pValue,breaks=c(-Inf, 0.05, Inf), label=c("+", "")) 
  
  list.frames.thermal[[i]] <- LM.run.pft
  df.pft.thermal <- rbind(LM.run.pft, df.pft.thermal)
}

# We used a relatively smaller threshold for the minimum sampling number for a bin pair (n >= 10), leading to some outliers in terms of gamma
# Gamma lie outside the interval formed by the 5 and 95 percentiles will be considered as potential outliers
df_fig2b <- df.pft.thermal %>% dplyr::select(reg.coef, Group) %>% 
  dplyr::rename(gammaT = reg.coef, PFT = Group)%>%
  dplyr::group_by(PFT) %>% 
  dplyr::filter(gammaT > quantile(gammaT, 0.05, na.rm = TRUE),
                gammaT < quantile(gammaT, 0.95, na.rm = TRUE))

temp <- LMM.run %>% dplyr::select(reg.coef) %>% 
  mutate(PFT = "ALL") %>% dplyr::rename(gammaT = reg.coef)

df_fig2b <- rbind(temp, df_fig2b)
df_fig2b$PFT <- as.factor(df_fig2b$PFT)
pft.lm <- lm(gammaT ~ PFT, df_fig2b)
pft.av <- aov(pft.lm)
tukey.test <- HSD.test(pft.av, trt = 'PFT')
tapply(df_fig2b$gammaT,df_fig2b$PFT,mean)

Tukey_label.T <- data.frame(PFT = c("ALL" ,"CRO", "DBF", 
                                    "EBF", "ENF", "GRA",
                                    "MF","WET"),
                            Letter = c("cd", "a", "b", 
                                       "d", "bc", "d", 
                                       "cd", "b"),
                            gammaT = c(0.4075306, 0.7726447, 0.5684652, 
                                       0.3782239, 0.5258282,0.3445662, 
                                       0.4089037, 0.5780809))
df_fig2b$PFT2 <- factor(df_fig2b$PFT,
                        levels = c('CRO','WET','DBF',
                                   'ENF','MF','EBF','GRA','ALL'),ordered = TRUE)

fig2b <- ggplot(df_fig2b, aes(x= PFT2, y=gammaT, color = PFT2)) + 
  geom_boxplot(size = 1, outlier.color= NA) + 
  geom_jitter(position=position_jitter(0.15), size = 0.4, alpha = 0.3) + 
  geom_text(data = Tukey_label.T, aes(x = PFT, label = Letter), vjust=-12, size = 5, color = "black") + 
  ylab(expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')'))+
  geom_hline(yintercept = 0, color = "black", linetype = 2) + 
  scale_y_continuous(limit = c(-2, 4))+
  scale_color_npg()+
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank()) +
  guides(color = FALSE)

# 3 - Analysis: Partial correlation on single sites--------------------------------
SiteInfo <- read.csv(paste0(InputFolder, "SiteInfo.csv")) %>% drop_na(nYear) %>% 
  filter(!(PFT %in% c("SNOW")))
# select sites with observations more than 6 yrs during MODIS periods 
# 77 sites were left;"US-Cop" has no data after filtering
df_LongObs <- SiteInfo %>% filter(nYear >= 6 & nYearMODIS > 6, 
                                 !(PFT %in% c("WSH", "SH")), !(Site %in% c("US-Cop")))

df_SingleSite <- data.frame()

for (NameSite in df_LongObs$Site){
  
  site.clean <- df.main %>% filter(Site == NameSite) %>% 
    filter(fAPAR >= 0.3, LagT > 0, daytime_Tair > 0, LagT > 0, 
           LagAlpha >= 0.7, daytime_VPD < 2) %>% drop_na(Amax2000)
  
  pcorr_LagTair <- pcor.test(site.clean$Amax2000, 
                             site.clean$LagT, 
                             c(site.clean$LagPAR,
                               site.clean$daytime_Tair,
                               site.clean$fAPAR))[,1]
  pValue_LagTair <- pcor.test(site.clean$Amax2000, 
                              site.clean$LagT, 
                              c(site.clean$LagPAR,
                                site.clean$daytime_Tair,
                                site.clean$fAPAR))[,2]
  nSample_LagTair <- pcor.test(site.clean$Amax2000, 
                               site.clean$LagT, 
                               c(site.clean$LagPAR,
                                 site.clean$daytime_Tair,
                                 site.clean$fAPAR))[,4]
  SiteRecord <- data.frame(Site = NameSite, pcorr = pcorr_LagTair, pValue = pValue_LagTair, 
                           PFT = site.clean$PET[1], nSample = nSample_LagTair)
  df_SingleSite <- rbind(SiteRecord, df_SingleSite)
}
# 92.20779% is positive (71 out of 77)
# nrow(df_SingleSite[df_SingleSite$pcorr>0,])/nrow(df_SingleSite)
df_SingleSite <- df_SingleSite %>% 
  dplyr::rename(Group = PFT) %>%  
  mutate(PFT3 = case_when(Group == 12 ~ "CRO", Group == 4 ~ "DBF",
                         Group == 2 ~ "EBF", Group == 1 ~ "ENF",
                         Group == 5 ~ "MF", Group == 10 ~ "GRA",
                         Group == 11 ~ "WET"))
# color consistent with Fig2b
df_SingleSite$PFT2 <- factor(df_SingleSite$PFT,
                        levels = c('CRO','WET','DBF','ENF','MF','EBF','GRA','ALL'),ordered = TRUE)

fig2c <- ggplot(df_SingleSite, aes(x=reorder(Site, pcorr), y=pcorr, fill=PFT2)) +
  geom_bar(stat="identity", alpha=0.8,  width = 0.6)+
  xlab("Site") + 
  ylab(expression(paste("Partial ",italic("r")))) + 
  geom_hline(yintercept = 0, linetype="longdash", size = 0.8)+
  scale_fill_npg()+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(),
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(1, 'cm'),
        legend.title = element_blank(),
        legend.direction="horizontal",
        legend.position=c(0.2, 0.85),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# 4 - Figure Ext 1--------------------------------
# One example for the use of LMM under specific bins of Tair and fPAR
df.exa <- df.clean.cor %>% filter(bin_fAPAR == "(0.56,0.58]",bin_T == "(15,16]")

lmm.model <- lmer(FLUXNET~LagTair+(1|Site),data=df.exa)
# summary(lmm.model)
# sqrt(rsquared.glmm(lmm.model)$Marginal)
# sqrt(rsquared.glmm(lmm.model)$Conditional)

FigExt1 <- ggplot(df.exa, aes(x = LagTair, y = FLUXNET, colour=Site)) +
  labs(y= expression("A"["max,2000"]~'('*mu*'mol'~'m'^-2~'s'^-1*')'))+
  xlab(expression(paste(bar(T[air]), " (\u00B0C)")))+
  geom_point(shape = 16, size=3) +
  geom_abline(aes(intercept=`(Intercept)`, slope=LagTair), as.data.frame(t(fixef(lmm.model))))+
  scale_y_continuous(limits = c(-5,60))+
  scale_x_continuous(limits = c(0,25))+
  annotate("text", x=3, y=50, label= "Slope = 0.88, p-value < 0.001\nConditional r = 0.34\nMarginal r = 0.88", 
           hjust = 0) +
  ggtitle(bquote("fPAR = (0.56,0.58],"~T[air]~"= (15,16]"))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())

# 5 - Figure Ext 4--------------------------------

df_map <- merge(df_SingleSite, df_LongObs, by = "Site") %>% 
  dplyr::select(Site, Latitude, Longitude, pcorr, PFT, nYearMODIS, MeanTg, SdTg) %>% 
  dplyr::rename(PFT = PFT.y)


world_map <- borders("world", colour="black", fill="white")
df_map_ <- df_map %>% mutate(PFT = ifelse(PFT == "DBF" | PFT == "ENF" | PFT == "MF", "Forest", PFT))

FigExt4a <- ggplot() + 
  world_map +
  geom_point(data = df_map_, aes(x = Longitude, y = Latitude, fill = pcorr, shape = PFT), 
             size = 3, stroke = 0.5) +
  scale_fill_gradient2(low = "blue", high = "#DC0000B2", 
                       limits = c(-0.11001,0.8), midpoint = 0, oob=squish) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_y_continuous(breaks= c(-60, -30, 0, 30, 60, 90), limits=c(-60, 90)) +
  labs(tag = "a") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.title = element_blank()) +
  labs(x = "Longitude", y = "Latitude")


Fig4b <- ggplot(data = df_map_, aes(MeanTg, pcorr)) +
  geom_point(aes(shape = PFT), size = 1.5, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black")+
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  ylab("Partial correlation coefficient") +
  xlab(expression(paste("Average of ", bar(T[air]), " (\u00B0C)")))+
  scale_x_continuous(limits = c(5, 30), expand = c(0,0), breaks = seq(5, 30, 5)) +
  stat_cor(aes(label = paste(..p.label.., sep = "~`,`~")), label.x = 22, label.y = 0.7) +
  labs(tag = "b") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank())

Fig4c <- ggplot(data = df_map_, aes(SdTg, pcorr)) +
  geom_point(aes(shape = PFT), size = 1.5, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black")+
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  ylab("Partial correlation coefficient") +
  xlab(expression(paste("Standard deviation of ", bar(T[air]), " (\u00B0C)")))+
  scale_x_continuous(limits = c(0, 7), expand = c(0,0), breaks = seq(0, 7, 1)) +
  stat_cor(aes(label = paste(..p.label.., sep = "~`,`~")), label.x = 1, label.y = 0.7) +
  labs(tag = "c") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank())
