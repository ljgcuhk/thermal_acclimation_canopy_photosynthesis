rm(list=ls())
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, plyr, scales, gtools, lme4, lmerTest, ggsci, rcarbon, ppcor,
               psych, agricolae, cowplot, ggpubr)
InputFolder <- ".../0_file/"
df.main <- read.csv(paste0(InputFolder, "main.csv"))

df.clean <- df.main %>% filter(fAPAR >= 0.3, daytime_Tair > 0, LagT > 0, LagAlpha >= 0.7, 
                               daytime_VPD < 2,!(PET %in% c(99,12,7,9)))%>% drop_na(Amax) 
# Growth PPFD vs Growth Tair
FigS2a <- ggplot(df.clean, aes(LagT, LagPAR)) +
  geom_hex() +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 1500)) +
  stat_smooth(method = "lm", se = FALSE,
              color = "black", size = 1.3) +
  xlab(expression(paste(bar(T[air]), " (°C)")))+
  ylab(bquote(bar(PPFD)~"("*mu*"mol photons"~m^-2~s^-1*")")) +
  labs(tag = "a") +
  stat_regline_equation(label.x = 2, label.y = 620) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., 
                             sep = "~`,`~")), label.x = 2, label.y = 580) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Growth Tair vs Tair
FigS2b <- ggplot(df.clean, aes(daytime_Tair, LagT)) +
  geom_hex() +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       limits = c(0, 2000), oob = squish) +
  stat_smooth(method = "lm", se = FALSE,
              color = "black", size = 1.3) +
  ylab(expression(paste(bar(T[air]), " (°C)")))+
  xlab(expression(paste(T[air], " (°C)"))) +
  labs(tag = "b") +
  stat_regline_equation(label.x = 2, label.y = 33) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 2, label.y = 30) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Amax2000 vs VPD
FigS2c <- ggplot(df.clean, aes(daytime_VPD, Amax2000)) +
  geom_hex() +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       limits = c(0, 1000), oob = squish) +
  stat_smooth(method = "lm", se = FALSE,
              color = "black", size = 1.3) +
  xlab("VPD (kPa)")+
  ylab(bquote(A["max,2000"]~"("*mu*"mol"~m^-2~s^-1*")")) +
  labs(tag = "c") +
  stat_regline_equation(label.x = 0.8, label.y = 110) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0.8, label.y = 100) +
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
