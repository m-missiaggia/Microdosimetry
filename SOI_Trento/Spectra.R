library(stringr)
library(gridExtra)
library(maptools)
library(ggplot2)
library(ggalt)
library(ggthemes)
library(tibble)
library(viridis)
library(knitr)
library(plotly)
library(egg)
library(ggpubr)
library(GillespieSSA)
library(MASS)
library(rgl)
library(zoo)
library(Sim.DiffProc)
library(doParallel)
library(ggsci)
library(deSolve)
library(pracma)
library(diffeqr)
library(parallel)
library(tidymodels)
library(data.table)

# source("G:\\Other computers\\Il mio laptop\\Francesco\\Università\\Articoli\\Marta\\GSM2\\Dose Rate\\Code\\utilities_doserate.R")

theme_set(theme_bw()+theme(plot.title = element_text(size=20, color="black"),
                           axis.title.x = element_text(size=20, color="black"),
                           axis.title.y = element_text(size=20, color="black"),
                           axis.text.x = element_text(size=20, color="black"),
                           axis.text.y = element_text(size=20, color="black"),
                           legend.title = element_blank(),
                           legend.text = element_text(size=20, color="black")))

cb_a <- c("#0072B2","#000000","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_b <- c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")
cb_c <- c( "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

cb_nob <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

setwd("C:\\Users\\Utente\\Documents\\SOI_Trento\\14March2024\\")

df7 <- fread("70MeV_muscle_15nA_1.Spe")
count7 <- as.numeric(df7$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count7$channel <- 1:nrow(count7)
colnames(count7)[1] <- "count"
count7$count2 <- count7$count*count7$channel^2
count7$Type <- "15nA - 1 - 1403"
count7$count2 <- (count7$count2*max(count$count2[count7$channel > cut]))/max(count7$count2[count7$channel > cut])

cut_noise <- 21
count7$count[which(count7$channel < 27)] <- 0

m_MG <- 0.46
m_LG <- 0.444



bin<-10^(seq(log10(0.1),log10(1000),length.out=100))


compute_calibration_adc <- function(m, df_count, bin){

  y <- m*c(1:nrow(df_count))
  
  df_y<-data.frame(y = y, counts = df_count$count)
  
  count_h<-c()
  y_h<-c()
  for (b in 1:(length(bin)-1)) {
    count_h[b] <- sum(df_y$counts[(df_y$y > bin[b]) & (df_y$y <= bin[b+1])])
    y_h[b]<-(bin[b]+bin[b+1])/2
  }
  
  fy_h <- data.frame(y = y_h, fy = count_h)
  
  return(fy_h)
}

compute_spectra <- function(fy_h, bin){
  
  hist <- fy_h %>% 
    mutate(BinWidth = diff(bin))
  
  hist$fy_bw <- hist$fy/hist$BinWidth
  
  B <- 1/diff(log10(hist$BinWidth))[1]
  
  # B<-1/0.05
  C <- log(10)*diff(log10(hist$BinWidth))[1]
  
  hist$fy_bw_norm <- hist$fy_bw/(C*sum(hist$y*hist$fy_bw))
  # hist$fy_norm <- (hist$counts)/(C*sum(hist$y*hist$count))
  
  hist$yfy <- hist$fy_bw_norm*hist$y
  
  hist$yfy_norm <- hist$yfy/(C*sum(hist$y*hist$yfy))
  
  hist$ydy <- hist$yfy_norm*hist$y
  
  hist_M <- fy_h %>% 
    mutate(BinWidth = diff(bin))
  
  hist_M$fy_bw<-hist_M$fy/hist_M$BinWidth
  hist_M$fy_bw_norm<-hist_M$fy_bw/(C*sum(hist_M$y*hist_M$fy))
  
  hist_M$yfy<-hist_M$fy_bw*hist_M$y
  hist_M$yfy_norm<-hist_M$yfy/(C*sum(hist_M$y*hist_M$yfy))
  hist_M$ydy<-hist_M$fy_bw*hist_M$y*hist_M$y
  
  #yF value
  yF<-sum(hist_M$BinWidth*hist_M$yfy)/sum(hist_M$BinWidth*hist_M$fy_bw)
  
  #yD value
  yD<-sum(hist_M$BinWidth*hist_M$ydy)/sum(hist_M$BinWidth*hist_M$yfy)
  
  y0<-(150)^2
  
  ystar<-(y0*sum(((1-exp(-(hist$y^2)/(y0)))*hist$BinWidth*hist$fy_bw)))/(yF*sum(hist$BinWidth*hist$fy_bw))
  
  return(list(histogram = hist, yF = yF, yD = yD, ystar = ystar))
}


fy_h <- compute_calibration_adc(m_MG, count7, bin)

df_plot <- data.frame(y = fy_h$y, fy = fy_h$fy)

p_adc <- df_plot %>% 
  ggplot(aes(x=y,y=fy)) +
  geom_step(linewidth = 1) +
  scale_x_log10(label = scientific_10) +
  xlab("y [keV/um]") +
  ylab("f(y)")
ggplotly(p_adc)

spectra <- compute_spectra(fy_h, bin)
  
p_spec <- spectra$histogram %>% 
  ggplot(aes(y,ydy)) + 
  geom_step(linewidth = 1) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) + 
  ylab("yd(y)")
ggplotly(p_spec)






