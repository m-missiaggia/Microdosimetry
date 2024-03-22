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

compute_calibration_adc <- function(m, df_count, bin, mv){
  
  y <- m*df_count$mv
  
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
  
  ySK<-sum(hist_M$BinWidth*hist_M$ydy*hist_M$y)/sum(hist_M$BinWidth*hist_M$yfy)
  
  yK<-sum(hist_M$BinWidth*hist_M$ydy*hist_M$y*hist_M$y)/sum(hist_M$BinWidth*hist_M$yfy)
  
  y0<-(150)^2
  
  ystar<-(y0*sum(((1-exp(-(hist$y^2)/(y0)))*hist$BinWidth*hist$fy_bw)))/(yF*sum(hist$BinWidth*hist$fy_bw))
  
  return(list(histogram = hist, yF = yF, yD = yD, ystar = ystar, ySK = ySK, yK = yK))
}

select_gain <- function(ion, gain){
  
  if(ion == "H"){
    if(gain == "M"){
      m <- 0.46
    }else if(gain == "L"){
      m <- 0.444
    }else if(gain == "C"){
      m <- 0.866
    }else{
      return(-1)
    }
  }else if(ion == "C"){
    if(gain == "M"){
      m <- 0.866
    }else if(gain == "L"){
      m <- 0.866
    }else{
      return(-1)
    }
  }
  
  return(m)
}

medium <- read.csv(file = "C:/Users/Utente/Documents/SOI_Trento/silicon_linearization_med_gain_H.csv")
low <- read.csv(file = "C:/Users/Utente/Documents/SOI_Trento/silicon_linearization_low_gain_H.csv")


ion <- "H"
cut <- 1
bin<-10^(seq(log10(1),log10(10000),length.out=100))
gain <- "M"
mv <- medium

m <- select_gain(ion, gain)

df_1 <- fread("C:/Users/Utente/Documents/SOI_Trento/18March2024/70MeV_3-8cmRW3_protonedge_50nA_medium.Spe")

count <- as.numeric(df_1$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame() %>% 
  mutate(mv = mv$INPUT..mV.)
count <- count %>% 
  mutate(channel = 1:nrow(count))
colnames(count)[1] <- "count"
count$count2 <- count$count*count$channel^2
# count$Type <- list_files[1]
count$Type <- "Muscle"

count$count[which(count$channel < cut)] <- 0

fy_h <- compute_calibration_adc(m, count, bin, mv)

df_plot <- data.frame(y = fy_h$y, fy = fy_h$fy)

spectra_s <- compute_spectra(fy_h, bin)$histogram %>% 
  mutate(Cumulative = cumsum(ydy)/sum(ydy)) %>% 
  mutate(Type = "Edge")

p_spec <- spectra_s %>% 
  ggplot(aes(y, ydy)) + 
  geom_step(linewidth = 1) +
  # scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #   labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_log10() +
  scale_color_manual(values = c25) +
  ylab("yd(y)")
ggplotly(p_spec)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

setwd("C:\\Users\\Utente\\Documents\\SOI_Trento\\18March2024/")

df <- fread("70MeV_3-8cmRW3_10nA_Prova1.Spe")
# noise <- fread("noise1.Spe")

sum(as.numeric(df$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(df$`$SPEC_ID:`[9]," ")[[1]][1]) 
# sum(as.numeric(noise$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(noise$`$SPEC_ID:`[9]," ")[[1]][1]) 

# sum(as.numeric(df$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(df$`$SPEC_ID:`[9]," ")[[1]][1]) - sum(as.numeric(noise$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(noise$`$SPEC_ID:`[9]," ")[[1]][1])


df <- fread("70MeV_3-5cmRW3_SOI.Spe")
noise <- fread("noise2.Spe")

sum(as.numeric(df$`$SPEC_ID:`[(12+36):4107]))/as.numeric(strsplit(df$`$SPEC_ID:`[9]," ")[[1]][1]) 
sum(as.numeric(noise$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(noise$`$SPEC_ID:`[9]," ")[[1]][1]) 

sum(as.numeric(df$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(df$`$SPEC_ID:`[9]," ")[[1]][1]) - sum(as.numeric(noise$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(noise$`$SPEC_ID:`[9]," ")[[1]][1])

df <- fread("70MeV_3-8cmRW3_protonedge_50nA_medium.Spe")
count <- as.numeric(df$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame() 
count <- count %>% 
  mutate(channel = 1:nrow(count))
colnames(count)[1] <- "count"
count$count2 <- count$count*count$channel^2
count$Type <- "10nA"

p <- count %>% ggplot(aes(channel,count2)) +
  geom_step(linewidth = 1) +
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks(sides = "b")
ggplotly(p)

p <- count %>% 
  filter(channel > 36) %>% 
  ggplot(aes(channel,count2)) +
  geom_step(linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(labels = label_comma()) +
  xlab("Channel") + 
  ylab("Counts") 
ggplotly(p)

df <- fread("70MeV_3-8cmRW3_10nA.Spe")
count <- as.numeric(df$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame() %>% 
  mutate(channel = 1:nrow(count))
colnames(count)[1] <- "count"
count$count2 <- count$count*count$channel^2
count$Type <- "10nA - 3.8"

df2 <- fread("70MeV_3-8cmRW3_50nA.Spe")
count2 <- as.numeric(df2$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count2$channel <- 1:nrow(count2)
colnames(count2)[1] <- "count"
count2$count2 <- count2$count*count2$channel^2
count2$Type <- "50nA - 3.8"
count2$count2 <- (count2$count2*max(count$count2[count2$channel > 36]))/max(count2$count2[count2$channel > 36])

sum(count2$count[which(count2$channel > 36)])/532

df3 <- fread("70MeV_3-9cmRW3_50nA.Spe")
count3 <- as.numeric(df3$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count3$channel <- 1:nrow(count3)
colnames(count3)[1] <- "count"
count3$count2 <- count3$count*count3$channel^2
count3$Type <- "50nA - 3.9"
count3$count2 <- (count3$count2*max(count$count2[count3$channel > 36]))/max(count3$count2[count3$channel > 36])

p <- count %>% rbind(count2, count3) %>% 
  filter(channel > 36) %>% 
  ggplot(aes(channel, count2, color = Type)) +
  geom_step(linewidth = 1) +
  scale_x_log10() +
  scale_color_manual(values = cb_b) +
  geom_vline(xintercept = 570, color = "maroon", linewidth = 1) +
  xlab("Channel") + 
  ylab("Counts")
ggplotly(p)
p





################à

setwd("C:\\Users\\Utente\\Documents\\SOI_Trento\\15March2024/")

df <- fread("70MeV_3-8cmRW3_lowgain.Spe")
count <- as.numeric(df$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count <- count %>% 
  mutate(channel = 1:nrow(count))
colnames(count)[1] <- "count"
count$count2 <- count$count*count$channel^2
count$Type <- "10nA new"

sum(count$count[which(count$channel > 20)])/484

p <- count %>% ggplot(aes(channel,count)) +
  geom_step(linewidth = 1) +
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks(sides = "b")
ggplotly(p)

p <- count %>% rbind(count2) %>% 
  ggplot(aes(channel, count, color = Type)) +
  geom_step(linewidth = 1) +
  scale_x_log10() + 
  annotation_logticks(sides = "b")
ggplotly(p)


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}





p <- count %>% rbind(count2) %>% 
  ggplot(aes(channel,count)) +
  geom_step(linewidth = 1) +
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks(sides = "b")
ggplotly(p)



as.numeric(df$`$SPEC_ID:`[12:4107]) %>% as.data.frame() %>% 
  ggplot(aes(.)) + geom_histogram() + scale_x_log10()

br<-10^seq(log10(10^0),log10(max(count$.)),length.out = 30)

hist_<-hist(count$.[which(count$. > bin_start & count$. < bin_end)],breaks=br) #use hist() and specify your breaks
hist<-data.frame(count=hist_$counts,x=0.5*(br[-1]+br[-length(br)]),xmin=br[-length(br)],xmax=br[-1])
hist$BinWidth<-abs(hist$xmax - hist$xmin)
hist$dose <- hist$x*hist$count

hist %>% ggplot(aes(x,count)) + geom_line(linewidth = 1) + scale_x_log10()
hist %>% ggplot(aes(x,dose)) + geom_line(linewidth = 1) + scale_x_log10()


