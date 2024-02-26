rm(list=ls())
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, zoo, cowplot)
InputFolder <- ".../0_file/"

df_ext <- read.csv(paste0(InputFolder, "time_scale/Amax2000_EBF_ext_TimeScale.csv")) %>% 
  arrange(ts) %>% mutate(mean.r = sqrt(mean.r2), sd.r = sqrt(sd.r2))
df_evi <- read.csv(paste0(InputFolder, "time_scale/EVI_EBF_TimeScale.csv")) %>% 
  arrange(ts) %>% mutate(mean.r = sqrt(mean.r2), sd.r = sqrt(sd.r2))

window_size <- 5

df_ext$mean.r_smooth <- ifelse(is.na(rollapply(df_ext$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                               df_ext$mean.r, 
                               rollapply(df_ext$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))
df_evi$mean.r_smooth <- ifelse(is.na(rollapply(df_evi$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                               df_evi$mean.r, 
                               rollapply(df_evi$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))

T_max <- max(df_evi$mean.r_smooth)

FigS1a <- ggplot(df_ext, aes(ts, mean.r_smooth)) + 
  geom_line(size = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  scale_x_continuous(limits = c(0,181), expand = c(0,0), breaks = c(0,30,60,90,120,150,180))+
  ggtitle(bquote('A'[max*',2000']*','~tau~'= ' ~ 'N/A')) + 
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

FigS1b <- ggplot(df_evi, aes(ts, mean.r_smooth)) + 
  geom_line(size = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  ggtitle(bquote('EVI, '~tau~'= ' ~ '13 d')) + 
  geom_hline(yintercept = T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = 13, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())
