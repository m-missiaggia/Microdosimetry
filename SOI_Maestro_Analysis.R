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


setwd("C:\\Users\\Utente\\Documents\\SOI_Trento\\13March2024\\")

cut <- 22

df <- fread("70MeV_muscle_15nA.Spe")
count <- as.numeric(df$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count <- count %>% 
  mutate(channel = 1:nrow(count))
colnames(count)[1] <- "count"
count$count2 <- count$count*count$channel^2
count$Type <- "15nA"

df2 <- fread("70MeV_muscle_40nA.Spe")
count2 <- as.numeric(df2$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count2$channel <- 1:nrow(count2)
colnames(count2)[1] <- "count"
count2$count2 <- count2$count*count2$channel^2
count2$Type <- "40nA"
count2$count2 <- (count2$count2*max(count$count2[count2$channel > cut]))/max(count2$count2[count2$channel > cut])

df3 <- fread("70MeV_muscle_30nA_1.Spe")
count3 <- as.numeric(df3$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count3$channel <- 1:nrow(count3)
colnames(count3)[1] <- "count"
count3$count2 <- count3$count*count3$channel^2
count3$Type <- "30nA"
count3$count2 <- (count3$count2*max(count$count2[count3$channel > cut]))/max(count3$count2[count3$channel > cut])

df4 <- fread("70MeV_muscle_30nA_2.Spe")
count4 <- as.numeric(df4$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count4$channel <- 1:nrow(count4)
colnames(count4)[1] <- "count"
count4$count2 <- count4$count*count4$channel^2
count4$Type <- "30nA - 2"
count4$count2 <- (count4$count2*max(count$count2[count4$channel > cut]))/max(count4$count2[count4$channel > cut])

df5 <- fread("70MeV_muscle_30nA_3.Spe")
count5 <- as.numeric(df5$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count5$channel <- 1:nrow(count5)
colnames(count5)[1] <- "count"
count5$count2 <- count5$count*count5$channel^2
count5$Type <- "30nA - 3"
count5$count2 <- (count5$count2*max(count$count2[count5$channel > cut]))/max(count5$count2[count5$channel > cut])

df6 <- fread("70MeV_muscle_30nA_4.Spe")
count6 <- as.numeric(df6$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count6$channel <- 1:nrow(count6)
colnames(count6)[1] <- "count"
count6$count2 <- count6$count*count6$channel^2
count6$Type <- "30nA - 4"
count6$count2 <- (count6$count2*max(count$count2[count6$channel > cut]))/max(count6$count2[count6$channel > cut])

p <- count %>% rbind(count2, count3, count4, count5, count6) %>% 
  filter(channel > cut) %>% 
  ggplot(aes(channel, count2, color = Type)) +
  geom_step(linewidth = 1) +
  scale_x_log10() +
  scale_color_manual(values = cb_b) +
  # geom_vline(xintercept = 570, color = "maroon", linewidth = 1) +
  xlab("Channel") + 
  ylab("Counts")
ggplotly(p)
p



sum(as.numeric(df4$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(df4$`$SPEC_ID:`[9]," ")[[1]][1]) 
sum(as.numeric(df4$`$SPEC_ID:`[(12 + 21):4107]))/as.numeric(strsplit(df4$`$SPEC_ID:`[9]," ")[[1]][2]) 


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


