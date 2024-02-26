rm(list=ls())
#install.packages("pacman") implement it if "pacman" package is not installed
pacman::p_load(tidyverse, cowplot, zoo)
InputFolder <- ".../0_file/"

# 1 - read data----------------------------------------------------------------
df.main <- data.frame()
Metric <- "Amax2000"

for (PFT in c("CRO", "DBF", "ENF", "GRA", "WET", "ALL")){
  f <- list.files(paste0(InputPath, "time_scale/"), pattern= paste0(Metric, "_", PFT))
  df <- read.csv(paste0(InputPath, "time_scale/", f)) %>% arrange(ts) %>% 
    mutate(mean.r = sqrt(mean.r2), sd.r = sqrt(sd.r2), PFT = PFT, Metric = Metric)
  df.main <- rbind(df.main, df)
}

# 2 - time scale for PFTs and all sites (ALL)-----------------------------------
#CRO: Cropland
df.Amax.CRO <- df.main %>% filter(PFT == "CRO")
# Calculate the moving average
df.Amax.CRO$mean.r_smooth <- ifelse(is.na(rollapply(df.Amax.CRO$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                                    df.Amax.CRO$mean.r, 
                                    rollapply(df.Amax.CRO$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))
CRO_T_max <- max(df.Amax.CRO$mean.r_smooth)
CRO_T_opt <- df.Amax.CRO$ts[df.Amax.CRO$mean.r_smooth == CRO_T_max]

Fig1_a <- ggplot(df.Amax.CRO, aes(ts, mean.r_smooth)) + 
  geom_line(linewidth = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  ggtitle(bquote('CRO, '~tau~'= ' ~ .(CRO_T_opt) ~ ' d')) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  geom_hline(yintercept = CRO_T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = CRO_T_opt, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

# DBF: deciduous broadleaf forests
df.Amax.DBF <- df.main %>% filter(PFT == "DBF")
# Calculate the moving average
df.Amax.DBF$mean.r_smooth <- ifelse(is.na(rollapply(df.Amax.DBF$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                                    df.Amax.DBF$mean.r, 
                                    rollapply(df.Amax.DBF$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))

DBF_T_max <- max(df.Amax.DBF$mean.r_smooth)
DBF_T_opt <- df.Amax.DBF$ts[df.Amax.DBF$mean.r_smooth == DBF_T_max]

Fig1_b <- ggplot(df.Amax.DBF, aes(ts, mean.r_smooth)) + 
  geom_line(linewidth = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  ggtitle(bquote('DBF, '~tau~'= ' ~ .(DBF_T_opt) ~ ' d')) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  geom_hline(yintercept = DBF_T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = DBF_T_opt, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

# ENF: evergreen needleleaf forests
df.Amax.ENF <- df.main %>% filter(PFT == "ENF")
# Calculate the moving average
df.Amax.ENF$mean.r_smooth <- ifelse(is.na(rollapply(df.Amax.ENF$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                                    df.Amax.ENF$mean.r, 
                                    rollapply(df.Amax.ENF$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))

ENF_T_max <- max(df.Amax.ENF$mean.r_smooth)
ENF_T_opt <- df.Amax.ENF$ts[df.Amax.ENF$mean.r_smooth == ENF_T_max]

Fig1_c <- ggplot(df.Amax.ENF, aes(ts, mean.r_smooth)) + 
  geom_line(linewidth = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  ggtitle(bquote('ENF, '~tau~'= ' ~ .(ENF_T_opt) ~ ' d')) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  geom_hline(yintercept = ENF_T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = ENF_T_opt, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

# GRA: grasslands
df.Amax.GRA <- df.main %>% filter(PFT == "GRA")
# Calculate the moving average
df.Amax.GRA$mean.r_smooth <- ifelse(is.na(rollapply(df.Amax.GRA$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                                    df.Amax.GRA$mean.r, 
                                    rollapply(df.Amax.GRA$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))

GRA_T_max <- max(df.Amax.GRA$mean.r_smooth)
GRA_T_opt <- df.Amax.GRA$ts[df.Amax.GRA$mean.r_smooth == GRA_T_max]

Fig1_d <- ggplot(df.Amax.GRA, aes(ts, mean.r_smooth)) + 
  geom_line(linewidth = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  ggtitle(bquote('GRA, '~tau~'= ' ~ .(GRA_T_opt) ~ ' d')) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  geom_hline(yintercept = GRA_T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = GRA_T_opt, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

# WET: wetlands
df.Amax.WET <- df.main %>% filter(PFT == "WET")
# Calculate the moving average
df.Amax.WET$mean.r_smooth <- ifelse(is.na(rollapply(df.Amax.WET$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                                    df.Amax.WET$mean.r, 
                                    rollapply(df.Amax.WET$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))

WET_T_max <- max(df.Amax.WET$mean.r_smooth)
WET_T_opt <- df.Amax.WET$ts[df.Amax.WET$mean.r_smooth == WET_T_max]

Fig1_e <- ggplot(df.Amax.WET, aes(ts, mean.r_smooth)) + 
  geom_line(linewidth = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  ggtitle(bquote('WET, '~tau~'= ' ~ .(WET_T_opt) ~ ' d')) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  geom_hline(yintercept = WET_T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = WET_T_opt, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

# All: All sites
df.Amax.ALL <- df.main %>% filter(PFT == "ALL")
# Calculate the moving average
df.Amax.ALL$mean.r_smooth <- ifelse(is.na(rollapply(df.Amax.ALL$mean.r, width = window_size, FUN = mean, align = "center", fill = NA)), 
                                    df.Amax.ALL$mean.r, 
                                    rollapply(df.Amax.ALL$mean.r, width = window_size, FUN = mean, align = "center", fill = NA))

ALL_T_max <- max(df.Amax.ALL$mean.r_smooth)
ALL_T_opt <- df.Amax.ALL$ts[df.Amax.ALL$mean.r_smooth == ALL_T_max]

Fig1_f <- ggplot(df.Amax.ALL, aes(ts, mean.r_smooth)) + 
  geom_line(linewidth = 1, color = "black") + 
  xlab("Number of days") + 
  ylab(expression(italic(r))) + 
  ggtitle(bquote('ALL, '~tau~'= ' ~ .(ALL_T_opt) ~ ' d')) + 
  scale_x_continuous(limits = c(0,61), expand = c(0,0), breaks = c(0,10,20,30,40,50,60))+
  geom_hline(yintercept = ALL_T_max, linetype="dashed", color = "black")+
  geom_vline(xintercept = ALL_T_opt, linetype="dashed", color = "black")+
  theme_bw()+ 
  theme(plot.title = element_text(size = 11,hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank())
