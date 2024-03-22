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

setwd("C:\\Users\\Utente\\Documents\\SOI_Trento\\14March2024\\")

df7 <- fread("70MeV_muscle_15nA_1.Spe")
count7 <- as.numeric(df7$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count7$channel <- 1:nrow(count7)
colnames(count7)[1] <- "count"
count7$count2 <- count7$count*count7$channel^2
count7$Type <- "15nA - 1 - 1403"
count7$count2 <- (count7$count2*max(count$count2[count7$channel > cut]))/max(count7$count2[count7$channel > cut])

df8 <- fread("70MeV_muscle_15nA_2.Spe")
count8 <- as.numeric(df8$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count8$channel <- 1:nrow(count8)
colnames(count8)[1] <- "count"
count8$count2 <- count8$count*count8$channel^2
count8$Type <- "15nA - 2 - 1403"
count8$count2 <- (count8$count2*max(count$count2[count8$channel > cut]))/max(count8$count2[count8$channel > cut])

df9 <- fread("70MeV_muscle_15nA_3.Spe")
count9 <- as.numeric(df9$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count9$channel <- 1:nrow(count9)
colnames(count9)[1] <- "count"
count9$count2 <- count9$count*count9$channel^2
count9$Type <- "15nA - 3 - 1403"
count9$count2 <- (count9$count2*max(count$count2[count9$channel > cut]))/max(count9$count2[count9$channel > cut])

df10 <- fread("70MeV_muscle_15nA_4.Spe")
count10 <- as.numeric(df10$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count10$channel <- 1:nrow(count10)
colnames(count10)[1] <- "count"
count10$count2 <- count10$count*count10$channel^2
count10$Type <- "15nA - 4 - 1403"
count10$count2 <- (count10$count2*max(count$count2[count10$channel > cut]))/max(count10$count2[count10$channel > cut])

df11 <- fread("70MeV_muscle_15nA_5.Spe")
count11 <- as.numeric(df11$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count11$channel <- 1:nrow(count11)
colnames(count11)[1] <- "count"
count11$count2 <- count11$count*count11$channel^2
count11$Type <- "15nA - 5 - 1403"
count11$count2 <- (count11$count2*max(count$count2[count11$channel > cut]))/max(count11$count2[count11$channel > cut])

df12 <- fread("70MeV_muscle_15nA_6.Spe")
count12 <- as.numeric(df12$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count12$channel <- 1:nrow(count12)
colnames(count12)[1] <- "count"
count12$count2 <- count12$count*count12$channel^2
count12$Type <- "15nA - 6 - 1403"
count12$count2 <- (count12$count2*max(count$count2[count12$channel > cut]))/max(count12$count2[count12$channel > cut])

df13 <- fread("70MeV_muscle_15nA_7.Spe")
count13 <- as.numeric(df13$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count13$channel <- 1:nrow(count13)
colnames(count13)[1] <- "count"
count13$count2 <- count13$count*count13$channel^2
count13$Type <- "15nA - 7 - 1403"
count13$count2 <- (count13$count2*max(count$count2[count13$channel > cut]))/max(count13$count2[count13$channel > cut])

df14 <- fread("70MeV_muscle_15nA_8.Spe")
count14 <- as.numeric(df14$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count14$channel <- 1:nrow(count14)
colnames(count14)[1] <- "count"
count14$count2 <- count14$count*count14$channel^2
count14$Type <- "15nA - 8 - 1403"
count14$count2 <- (count14$count2*max(count$count2[count14$channel > cut]))/max(count14$count2[count14$channel > cut])


p <- count %>% rbind(count7, count8, count9, count10, count11) %>% 
  filter(channel > cut) %>% 
  ggplot(aes(channel, count2, color = Type)) +
  geom_step(linewidth = 1) +
  scale_x_log10() +
  scale_color_manual(values = cb_b) +
  # geom_vline(xintercept = 570, color = "maroon", linewidth = 1) +
  xlab("Channel") + 
  ylab("Counts") +
  geom_vline(xintercept = 570, color = "maroon", linewidth = 1)
ggplotly(p)
p

count_sum <- as.numeric(df9$`$SPEC_ID:`[12:4107]) + as.numeric(df10$`$SPEC_ID:`[12:4107]) +
  as.numeric(df11$`$SPEC_ID:`[12:4107]) + as.numeric(df12$`$SPEC_ID:`[12:4107]) + 
  as.numeric(df13$`$SPEC_ID:`[12:4107]) +   as.numeric(df14$`$SPEC_ID:`[12:4107])  %>% 
  as.data.frame()
count_sum$channel <- 1:nrow(count_sum)
colnames(count_sum)[1] <- "count"
count_sum$count2 <- count_sum$count*count_sum$channel^2
count_sum$Type <- "15nA - sum"
count_sum$count2 <- (count_sum$count2*max(count$count2[count_sum$channel > cut]))/max(count_sum$count2[count_sum$channel > cut])

psum <- count %>% rbind(count_sum) %>% 
  filter(channel > cut) %>% 
  ggplot(aes(channel, count2, color = Type)) +
  geom_step(linewidth = 1) +
  scale_x_log10() +
  scale_color_manual(values = cb_b) +
  # geom_vline(xintercept = 570, color = "maroon", linewidth = 1) +
  xlab("Channel") + 
  ylab("Counts") +
  geom_vline(xintercept = 570, color = "maroon", linewidth = 1)
ggplotly(psum)

count_sum$count[which(count_sum$channel > cut)] %>% sum()




############1603

setwd("C:\\Users\\Utente\\Documents\\SOI_Trento\\16March2024\\")

cut <- 16

df <- fread("70MeV_muscle_15nA_.Spe")
count <- as.numeric(df$`$SPEC_ID:`[12:4107]) %>% 
  as.data.frame()
count <- count %>% 
  mutate(channel = 1:nrow(count))
colnames(count)[1] <- "count"
count$count2 <- count$count*count$channel^2
count$Type <- "15nA - muscle"

sum(as.numeric(df$`$SPEC_ID:`[12:4107]))/as.numeric(strsplit(df$`$SPEC_ID:`[9]," ")[[1]][1]) 
sum(as.numeric(df$`$SPEC_ID:`[(12 + 21):4107]))/as.numeric(strsplit(df$`$SPEC_ID:`[9]," ")[[1]][2]) 

p <- count %>% 
  filter(channel > cut) %>% 
  ggplot(aes(channel, count2, color = Type)) +
  geom_step(linewidth = 1) +
  scale_x_log10() +
  scale_color_manual(values = cb_b) +
  # geom_vline(xintercept = 570, color = "maroon", linewidth = 1) +
  xlab("Channel") + 
  ylab("Counts")
ggplotly(p)






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


