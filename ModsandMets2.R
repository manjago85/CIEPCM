library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ncdf4)
library(raster)
#library(HistogramTools)
library(ggrepel)
library(transport)
library(fda)
library(fda.usc)
library(stringr)
library(tidyr)
library(doBy)
library(topsis)
library(PROMETHEE)
library(cowplot)
library(RColorBrewer)
library(viridis)
# library(sf)
# library(rgeos)

#-------------------------------------------------------
#FUNCTIONS
#-------------------------------------------------------

source("Pfunctions.R")

#-------------------------------------------------------
#Color palette
#-------------------------------------------------------
color.pal <- brewer.pal(12,"Paired")
#-------------------------------------------------------
#Model's data
#-------------------------------------------------------
f <- "model_info_2.csv"
mdsdf <- read.csv(f,
                  header = TRUE,
                  sep = ",")

#Area of interest
area.lmts <- extent(265,285, 4, 20)
area.tmpl <- raster(ext = area.lmts,res = c(2.5,2.5))
#Teleconnection's areas
nino1.2.lmts <- extent(270,280, -10, 0)
nino1.2.tmpl <- raster(ext = nino1.2.lmts,res = c(2.0,2.0)) 
nino3.lmts <- extent(210,270, -5, 5)
nino3.tmpl <- raster(ext = nino3.lmts,res = c(2.0,2.0))
nino3.4.lmts <- extent(190,240, -5, 5)
nino3.4.tmpl <- raster(ext = nino3.4.lmts,res = c(2.0,2.0))
nino4.lmts <- extent(160,210, -5, 5)
nino4.tmpl <- raster(ext = nino4.lmts,res = c(2.0,2.0))
atn.lmts <- extent(302.5,345, 5.5, 23.5)
atn.tmpl <- raster(ext = atn.lmts,res = c(2.0,2.0))
alltlc.lmts <- extent(160,345, -5, 5)
alltlc.tmpl <- raster(ext = alltlc.lmts,res = c(2.0,2.0))
tlc.lmts.list <- list(nino1.2.lmts,nino3.lmts,nino3.4.lmts,nino4.lmts,atn.lmts)
tlc.tmpl.list <- list(nino1.2.tmpl,nino3.tmpl,nino3.4.tmpl,nino4.tmpl,atn.tmpl)

#-----------------------------------------------------
#PARAMETERS
#-----------------------------------------------------
mdq <- 32 #Model quantity
tst <- 1:12 #Month indexes
prd <- length(tst) #Period
yrq <- 36 #Year quantity
#-----------------------------------------------------
#MATRICES
#-----------------------------------------------------
#List to store results
measures.names <- c("pr.mean","pr.sd",
                    "tas.mean","tas.sd",
                    "nino1.2","nino3","nino3.4","nino4",
                    "atn") 
measures.obs <- vector("list", 9)
names(measures.obs) <- measures.names 
model.list <- vector("list", mdq)
measures.mod <- vector("list", 9)
for (p in 1:9){
  measures.mod[[p]] <- model.list
}
names(measures.mod) <- measures.names 

#Vectors to store models names
pr.names <- character(mdq)
ta.names <- character(mdq)

#-----------------------------------------------------
#PRECIPITATION
#-----------------------------------------------------
d <- "/home/mgomez/Dropbox/Tesis/modelos_Almazroui/Precipitacion"
setwd(d)
#Reference data
obs <- process.raster("precip.mon.mean_v23.nc",
                      "precip",
                      "1979-01-01",
                      "2015-01-01",
                      area.lmts,
                      c(1,0),
                      area.tmpl)
#Index
ind <- month(getZ(obs))
#Mean
measures.obs$pr.mean <- group.raster(obs,ind,mean)@data@values
#pr.yearmean.obs <- raster::calc(obs,fun = mean,na.rm = TRUE)
#Standard deviation
measures.obs$pr.sd <- group.raster(obs,ind,sd)@data@values
#pr.yearsd.obs <- raster::calc(obs,fun = sd,na.rm = TRUE)
#----------------------------------------------
#Loop to generate deltas
for(k in 1:mdq){
print(k)
  ptt <- paste0("^",k,"-","\\w") #Defines name pattern
  f <- list.files(pattern = ptt) #Find file
  #mod.name <- str_extract(f,"\\_(.*)\\.")
  mod.name <- sub(".*pr_(.*)\\.nc$", "\\1", f)
  pr.names[k] <- mod.name
  #Process file
  mod <- process.raster(f,
                        "pr",
                        "1979-01-01",
                        "2015-01-01",
                        area.lmts,
                        c(86400,0),
                        area.tmpl)
  #Index
  ind <- month(getZ(mod))
  #-----------------------------------------
  #Mean
  mod.mean <- group.raster(mod,ind,mean)@data@values 
  #-----------------------------------------
  #Standard deviation
  mod.sd <- group.raster(mod,ind,sd)@data@values 

  #Store model's data
  measures.mod$pr.mean[[k]] <- mod.mean
  measures.mod$pr.sd[[k]] <- mod.sd
}
#-----------------------------------------------------
#TEMPERATURE
#-----------------------------------------------------
d <- "/home/mgomez/Dropbox/Tesis/modelos_Almazroui/Temperatura"

setwd(d)
#Reference data
obs <- process.raster("air.mon.mean.nc",
                      "air",
                      "1979-01-01",
                      "2015-01-01",
                      area.lmts,
                      c(1,0),
                      area.tmpl)
#Index
ind <- month(getZ(obs))
#Mean
measures.obs$tas.mean <- group.raster(obs,ind,mean)@data@values
#Standard deviation
measures.obs$tas.sd <- group.raster(obs,ind,sd)@data@values
#----------------------------------------------
#Loop to generate deltas
for(k in 1:mdq){
  print(k)
  ptt <- paste0("^",k,"-","\\w") #Defines name pattern
  f <- list.files(pattern = ptt) #Finds file
  #mod.name <- str_extract(f,"\\_(.*)\\_")
  mod.name <- sub(".*tmp_Amon_(.*)\\.nc$", "\\1", f)
  
  mod <- process.raster(f,
                        "tas",
                        "1979-01-01",
                        "2015-01-01",
                        area.lmts,
                        c(1,-273),
                        area.tmpl)
  #Index
  ind <- month(getZ(mod))
  #-----------------------------------------
  #Mean
  mod.mean <- group.raster(mod,ind,mean)@data@values 
  #-----------------------------------------
  #Standard deviation
  mod.sd <- group.raster(mod,ind,sd)@data@values  
  #Store model's data
  measures.mod$tas.mean[[k]] <- mod.mean
  measures.mod$tas.sd[[k]] <- mod.sd
}
#-----------------------------------------------------
#TELECONNECTIONS
#-----------------------------------------------------
d <- "/home/mgomez/Dropbox/Tesis/modelos_Almazroui/Precipitacion"
setwd(d)
#Reference data
pr.obs <- process.raster("precip.mon.mean_v23.nc",
                      "precip",
                      "1978-12-01",
                      "2014-12-01",
                      area.lmts,
                      c(1,0),
                      area.tmpl)
ind <- quarter(getZ(pr.obs),
               with_year = TRUE,
               fiscal_start = 12)
pr.tlc.obs <- group.raster(pr.obs,ind,mean)@data@values

for(t in 1:5){
print(t)
sst.obs <- process.raster("sst.mnmean_v5.nc",
                      "sst",
                      "1978-12-01",
                      "2014-12-01",
                      tlc.lmts.list[[t]],
                      c(1,0),
                      tlc.tmpl.list[[t]])
ind <- quarter(getZ(sst.obs),
               with_year = TRUE,
               fiscal_start = 12)
sst.tlc.obs <- group.raster(sst.obs,ind,mean)@data@values

#Time series
tlc.ts <- apply(sst.tlc.obs,2,mean,na.rm = TRUE)

p <- t + 4 #position in measure list
measures.obs[[p]] <- tlc.ptn(tlc.ts,
                              pr.tlc.obs,
                              sls = 48,
                              yrq = yrq)

#Loop to generate deltas
for(k in 1:mdq){
  print(k)
  ptt <- paste0("^",k,"-","\\w") #Defines name pattern
  f <- list.files(pattern = ptt) #Find file
  #mod.name <- str_extract(f,"\\_(.*)\\_")
  #pr.names[k] <- mod.name
  #Process file
  mod <- process.raster(f,
                        "pr",
                        "1978-12-01",
                        "2014-12-01",
                        area.lmts,
                        c(86400,0),
                        area.tmpl)
  #Index
  ind <- quarter(getZ(mod),
                 with_year = TRUE,
                 fiscal_start = 12)
  pr.tlc.mod <- group.raster(mod,ind,mean)@data@values
  #-----------------------------------------
  #Correlation
  tlc.mod <- tlc.ptn(tlc.ts,
                     pr.tlc.mod,
                     sls = 48,
                     yrq = yrq)
  #Store model's data
  measures.mod[[p]][[k]] <- tlc.mod
  }
}

#-----------------------------------------------------
#WASPAEF calculation and ranking
#-----------------------------------------------------

f <- "/home/mgomez/Dropbox/Tesis/"
setwd(f)
#save(measures.obs, file = "measuresobs2.Rdata")
#save(measures.mod, file = "measuresmod2.Rdata")
load("measuresobs2.Rdata")
load("measuresmod2.Rdata")

#Results list
measures.results <- vector("list", 10)
names(measures.results) <- c(measures.names,"enso")

for(p in 1:9){
  mea <- calc.measure(measures.obs[[p]],
                      measures.mod[[p]],
                      waspaef)
  measures.results[[p]] <- mea
}

#Combines ENSO results
measures.results[[10]] <- cbind(measures.results[[5]],
                                measures.results[[6]],
                                measures.results[[7]],
                                measures.results[[8]],
                                measures.results[[9]])

#Store 1-waspaef
measures.results.1 <- lapply(measures.results, 
                       FUN = function(l) 1-l)

#Calculation of deltas
measures.deltas <- lapply(measures.results, FUN = function(l){
apply(l,1,deltaf)  
}
)

deltasmat <- matrix(unlist(measures.deltas),nrow = mdq)
colnames(deltasmat) <- c(measures.names,"enso")
fmeasures.names <- c("pr.mean","pr.sd",
                     "tas.mean","tas.sd",
                     "enso","atn")
msq <- length(fmeasures.names)
deltasmat <- deltasmat[,fmeasures.names]
  

#Ranking - Euclidean distance
L2NM <- apply(deltasmat, 1,normvec) 
rank1 <- rank(L2NM)
#Ranking - TOPSIS
topsis.w <- rep(1,6)
topsis.s <- rep("-",6)
TPSS <- topsis(deltasmat, topsis.w, topsis.s)
rank2 <- TPSS[,3]
#Ranking - PROMETHEE
# PreferenceF
PreF<-as.data.frame(matrix(rep("Gaussian",mdq*msq),nrow = mdq))
colnames(PreF) <-  fmeasures.names
# PreferenceT
PreT<-as.data.frame(matrix(rep(0.2,mdq*msq),nrow = mdq))
colnames(PreT) <-  fmeasures.names
# IndifferenceT
IndT<-as.data.frame(matrix(rep(0.1,mdq*msq),nrow = mdq))
colnames(IndT) <-  fmeasures.names
#Weights
Weig<-as.data.frame(matrix(rep(1/6,mdq*msq),nrow = mdq))
colnames(Weig) <- fmeasures.names
# Min_Max
MiMa<-as.data.frame(matrix(rep("min",mdq*msq),nrow = mdq))
colnames(MiMa) <- fmeasures.names
#S_Gauss
gauss<- as.data.frame(matrix(rep(0,mdq*msq),nrow = mdq))
colnames(gauss) <- fmeasures.names

PMTH <-  PROMETHEE(as.data.frame(deltasmat),
                   PreF,
                   PreT,
                   IndT,
                   Weig,
                   MiMa,
                   gauss)
rank3 <- rank(-PMTH$PROMETHEE2)

qm <- 10
t1 <- which.minn(rank1,qm)
t2 <- which.minn(rank2,qm)
t3 <- which.minn(rank3,qm)
fmds <- intersect(intersect(t1,t2),t3)
fmds.names <- mdsdf %>%
  filter(nmod %in% fmds) %>%
  dplyr::select(ordmod_Identificador) %>%
  pull(ordmod_Identificador)

#Ranking Data frame
namesdf <- data.frame(name = mdsdf$Identificador,
                      variant = mdsdf$Variante) 

 fmeasures.names <- c("PRM","PRS","TPM","TPS","ENSO","TNA")
 colnames(deltasmat) <- fmeasures.names

rankdf <- apply(deltasmat,2,rank) %>%
  as.data.frame() %>% 
  mutate_all(as.numeric) %>% 
  bind_cols(namesdf) %>%
  relocate(name,variant) %>%
  bind_cols(L2N = rank1,TOPSIS = rank2,PROMETHEE = rank3) %>%
  mutate(nmod = factor(1:32)) %>% 
  mutate(ordmod = factor(1:32))

#Ranking chart
rankdf.l <- rankdf %>%
  filter(nmod %in% fmds) %>%
  mutate(identifier = str_c(ordmod,".",name,"\n",variant)) %>%
  mutate(identifier = factor(identifier,
                             levels = identifier)) %>% 
  pivot_longer(!c(identifier,nmod,ordmod,name,variant),
               names_to = "Measure",
               values_to = "Value") %>%
  mutate(Measure = factor(Measure,
                          levels = c(fmeasures.names,
                                     "L2N",
                                     "TOPSIS",
                                     "PROMETHEE")))

#Test
#library(ggplot2)
ggplot(rankdf.l, aes(x = Measure, 
                     y = identifier, 
                     fill = Value)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  geom_text(aes(label = Value),
            color = "white", 
            size = 3) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits=rev) +
  labs(x = NULL, 
       y = NULL,
       fill = "Ranking") +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin=margin(t=-10),
        axis.text.x = element_text(angle = 90,
                                   size = 6),
        axis.text.y = element_text(size = 6),
        plot.margin=unit(c(0,0,0,0),"mm")) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barwidth = 10,
                                barheight = 0.5)) +
  
  scale_fill_viridis(limits = c(1,32)) +
  
  coord_fixed()
ggsave("rank6.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)

#Boxplots
month.letters <- month.abb

#pr.mean
p1 <- bxp(measures.results.1[[1]],
          month.letters,
          fmds,
          lyt = c(2,3))
p1
ggsave("prmeanbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)
#pr.sd
p2 <- bxp(measures.results.1[[2]],
          month.letters,
          fmds,
          lyt = c(2,3))
p2
ggsave("prsdbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)
#tas.mean
p3 <- bxp(measures.results.1[[3]],
          month.letters,
          fmds,
          lyt = c(2,3))
p3
ggsave("tasmeanbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)
#tas.sd
p4 <- bxp(measures.results.1[[4]],
          month.letters,
          fmds,
          lyt = c(2,3))
p4
ggsave("tassdbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)

season.letters <- c("DJF","MAM","JJA","SON")
#season.letters <- c("DEF","MAM","JJA","SON")

#enso
p5 <- bxp(measures.results.1[[5]],
          season.letters,
          fmds,
          lyt = c(2,3))
p5
ggsave("enso12dbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)

p6 <- bxp(measures.results.1[[6]],
          season.letters,fmds,
          lyt = c(2,3))
p6
ggsave("enso3dbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)

p7 <- bxp(measures.results.1[[7]],
          season.letters,
          fmds,
          lyt = c(2,3))
p7
ggsave("enso34dbx12.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)

p8 <- bxp(measures.results.1[[8]],
          season.letters,
          fmds,
          lyt = c(2,3))
p8
ggsave("enso4dbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)

#atn
p9 <- bxp(measures.results.1[[9]],
          season.letters,
          fmds,
          lyt = c(2,3))
p9
ggsave("atndbx.png",
       dpi = "retina",
       units = "cm",
       width = 12,
       height = 8)
