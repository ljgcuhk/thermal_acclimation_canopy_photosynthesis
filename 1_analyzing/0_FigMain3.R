rm(list=ls())
#import functions used for statistics
source(paste0(InputFolder,"0_functions.R"))
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, scales, gtools, lme4, lmerTest, ggsci, rcarbon, ppcor,
               psych, agricolae, cowplot, ggpubr)

# 1 - read data----------------------------------------------------------------
InputFolder <- ".../0_file/"
df.main <- read.csv(paste0(InputFolder, "main.csv"))
df.site.info <- read.csv(paste0(InputFolder, "SiteInfo.csv"))
# the column "Model" indicates if the site is included in the analysis
# only the C3 and non-dry sites (71) with good estimation for Amax2000 (R2 > 0.5) are included 
model.site.list <- df.site.info %>% filter(Model == "Yes")
df.clean <- df.main %>% filter(Site %in% model.site.list$Site,
                               fAPAR >= 0.3, daytime_Tair > 0, LagT > 0,
                               LagAlpha >= 0.7, daytime_VPD < 2)

df.clean <-transform(df.clean, bin_fAPAR = cut(df.clean$fAPAR, c(seq(0.3, 1, by = 0.02), Inf), include.lowest = TRUE))
df.clean <-transform(df.clean, bin_T = cut(df.clean$daytime_Tair, c(seq(0, 35, by = 1), Inf), include.lowest = TRUE))
df.clean.Amax <- data.frame(FLUXNET = df.clean$Amax2000, 
                            LagTair = df.clean$LagT,
                            bin_fAPAR = df.clean$bin_fAPAR, 
                            bin_T = df.clean$bin_T,
                            Site = as.factor(df.clean$Site))

# v0 represents PFT-fiexed Vcmax_25C
# v1 represents LAI-scaled Vcmax_25C
# v2 represents optimality-based Vcmax_25C
df.clean.v0 <- data.frame(FLUXNET = df.clean$Amax2000_v0, 
                          LagTair = df.clean$LagT,
                          bin_fAPAR = df.clean$bin_fAPAR, 
                          bin_T = df.clean$bin_T,
                          Site = as.factor(df.clean$Site)) %>% 
                          drop_na(FLUXNET)

df.clean.v1 <- data.frame(FLUXNET = df.clean$Amax2000_v1, 
                          LagTair = df.clean$LagT,
                          bin_fAPAR = df.clean$bin_fAPAR, 
                          bin_T = df.clean$bin_T,
                          Site = as.factor(df.clean$Site)) %>% 
                          drop_na(FLUXNET)

df.clean.v2 <- data.frame(FLUXNET = df.clean$Amax2000_v2, 
                          LagTair = df.clean$LagT,
                          bin_fAPAR = df.clean$bin_fAPAR, 
                          bin_T = df.clean$bin_T,
                          Site = as.factor(df.clean$Site)) %>% 
                          drop_na(FLUXNET)

# 2 - Analysis: comparing obs and est-------------------------------------------
LMM.run.Amax <- data.frame(ddply(df.clean.Amax, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run.Amax$sig.flag <- cut(LMM.run.Amax$pValue,breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

LMM.run.v0 <- data.frame(ddply(df.clean.v0, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run.v0$sig.flag <- cut(LMM.run.v0$pValue,breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

LMM.run.v1 <- data.frame(ddply(df.clean.v1, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run.v1$sig.flag <- cut(LMM.run.v1$pValue,breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

LMM.run.v2 <- data.frame(ddply(df.clean.v2, .(bin_fAPAR,  bin_T), LMM_ALL))
LMM.run.v2$sig.flag <- cut(LMM.run.v2$pValue,breaks=c(-Inf, 0.05, Inf), label=c("\u2022", ""))

df0 <- LMM.run.Amax %>% dplyr::select(reg.coef, sig.flag, bin_fAPAR, bin_T) %>% 
  dplyr::rename(FLUXNET2015 = reg.coef, sig.flag.FLUXNET = sig.flag)
df1 <- LMM.run.v0 %>% dplyr::select(reg.coef, sig.flag, bin_fAPAR, bin_T) %>% 
  dplyr::rename(BESSv0 = reg.coef, sig.flag.v0 = sig.flag)
df2 <- LMM.run.v1 %>% dplyr::select(reg.coef, sig.flag, bin_fAPAR, bin_T) %>% 
  dplyr::rename(BESSv1 = reg.coef, sig.flag.v1 = sig.flag)
df3 <- LMM.run.v2 %>% dplyr::select(reg.coef, sig.flag, bin_fAPAR, bin_T) %>% 
  dplyr::rename(BESSv2 = reg.coef, sig.flag.v2 = sig.flag)

# we only compare the overlapping cells
df.merge <- merge(df0, df1, by = c("bin_fAPAR", "bin_T"), all.x = TRUE, all.y = FALSE)
df.merge <- merge(df.merge, df2, by = c("bin_fAPAR", "bin_T"), all.x = TRUE, all.y = FALSE)
df.merge <- merge(df.merge, df3, by = c("bin_fAPAR", "bin_T"), all.x = TRUE, all.y = FALSE) 
df.merge.final <- df.merge %>% drop_na()

df.his <- df.merge.final %>% 
  dplyr::select(FLUXNET2015, BESSv0, BESSv1, BESSv2) %>% 
  dplyr::rename(aFLUXNET2015 = FLUXNET2015, bBESS_PFT = BESSv0, 
                cBESS_LAI = BESSv1, dBESS_EEO = BESSv2) %>%
  drop_na() %>% 
  pivot_longer(cols = 1:4, names_to = "variable", values_to = "value") 

cols <- c("aFLUXNET2015" = "black","bBESS_PFT" = "#8491B4B2",
          "cBESS_LAI" = "#1F77B4", "dBESS_EEO" = "#FF7F0EFF")

stat_FLUXNET <- df.his %>% filter(variable == "aFLUXNET2015") %>% 
  dplyr::select(value) %>% summarise(median(value))
stat_BESS_PFT <- df.his %>% filter(variable == "bBESS_PFT") %>% 
  dplyr::select(value) %>% summarise(median(value))
stat_BESS_LAI <- df.his %>% filter(variable == "cBESS_LAI") %>% 
  dplyr::select(value) %>% summarise(median(value))
stat_BESS_EEO <- df.his %>% filter(variable == "dBESS_EEO") %>% 
  dplyr::select(value) %>% summarise(median(value))

Fig3_a <- ggplot(df.his, aes(color = variable)) + 
  geom_density(aes(value), size = 1 ,show.legend=FALSE) + 
  stat_density(aes(value),geom="line",position="identity", size = 0.8)+
  ylab("Probability density") + 
  xlab(expression(gamma[T]~'('*mu*'mol'~'m'^-2~'s'^-1~'\u00B0C'^-1*')'))+
  scale_color_manual(values = cols,
                     labels = c(expression("FLUXNET2015"), expression("BESS"[PFT]), 
                                expression("BESS"[LAI]), expression("BESS"[EEO])),
                     guide = guide_legend(label.position = "left", 
                                          label.hjust = 0))+
  geom_vline(aes(xintercept = stat_FLUXNET$`median(value)`),col='black',linetype = "dashed", size=0.8)+
  geom_vline(aes(xintercept = stat_BESS_PFT$`median(value)`),col="#8491B4B2",linetype = "dashed", size=0.8)+
  geom_vline(aes(xintercept = stat_BESS_LAI$`median(value)`),col="#1F77B4",linetype = "dashed", size=0.8)+
  geom_vline(aes(xintercept = stat_BESS_EEO$`median(value)`),col="#FF7F0EFF",linetype = "dashed", size=0.8)+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.direction="vertical",
    legend.position= c(0.2, 0.8),
    legend.background = element_blank(),
    legend.key = element_blank())

# 3 - Analysis: Kolmogorov–Smirnov test-------------------------------------------
ks_test_pft <- ks.test(df.merge.final$FLUXNET2015, df.merge.final$BESSv0)
ks_test_lai <- ks.test(df.merge.final$FLUXNET2015, df.merge.final$BESSv1)
ks_test_eeo <- ks.test(df.merge.final$FLUXNET2015, df.merge.final$BESSv2)

# Create a data frame for the results
ks_plot <- data.frame(
  comparison = c("aBESS_PFT", "bBESS_LAI", "cBESS_EEO"),
  ks_statistic = c(ks_test_pft$statistic, ks_test_lai$statistic, ks_test_eeo$statistic),
  p_value = c(ks_test_pft$p.value, ks_test_lai$p.value, ks_test_eeo$p.value),
  color = c("#8491B4B2", "#1F77B4", "#FF7F0EFF")
)

Fig3_b <- ggplot(ks_plot, aes(x = comparison, y = ks_statistic, fill = comparison)) +
  geom_bar(stat = "identity", width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = paste("p-value:", format.pval(p_value, digits = 2))),
            position = position_dodge(width = 0.6), vjust = -0.5) +
  scale_fill_manual(values = ks_plot$color) +
  labs(y = "K-S statistic") +
  scale_y_continuous(limits = c(0,0.7), expand = c(0, 0)) +
  scale_x_discrete(labels = c(expression("BESS"[PFT]), expression("BESS"[LAI]), expression("BESS"[EEO]))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank())
