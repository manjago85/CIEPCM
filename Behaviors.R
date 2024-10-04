library(dplyr)
library(ggplot2)
library(RandomFields)
library(fda.usc)
library(fda)
library(transport)
# library(proxy)
# library(refund)
# library(roahd)
# library(devtools)
#install_github("trevor-harris/kstat")
# library(kstat)
# library(verification)
# library(DescTools)
#library(hydroGOF)
library(HistogramTools)
# library(moments)
# library(refund)
# library(kSamples)
# library(roahd)

#-------------------------------------------------------
#Load functions
#-------------------------------------------------------

source("Pfunctions.R")

  #-------------------------------------------------------
  # Simulation scenarios
  #-------------------------------------------------------
  a <- c(-0.9,-0.3,0.3,0.9) #correlations
  b <- c(0.1,0.7,1,1.3,1.9) #standard deviation ratio
  c <- c(0,1,5) #bias
  d <- c(0,1) #Disturbed
  #Matrix to store combinations and method results
  Comb <- expand.grid(a,b,c,d);colnames(Comb) <- c("rho","lambda","epsilon","dbt")
  k <- nrow(Comb) #combinations
  #-----------------------------------------------------------
  #Define a model and spatio-temporal parameters
  #-----------------------------------------------------------
  #Creates desire model
  model <- RMmatern(nu = 1.5)
  #Defines spatial and temporal coordinates
  x <- 1:5
  y <- 1:5
  t <- 1:100
  xl <- length(x)
  yl <- length(y)
  tl <- length(t)
  #--------------------------------------------------------
  #Creates matrices to store models
  #--------------------------------------------------------
  mods.c <- matrix(data = NA,nrow = k ,ncol = xl*yl*10)
  mods.s <- matrix(data = NA,nrow = k ,ncol = xl*yl*tl)
  #-------------------------------------------------------
  #Creates spatio-temporal field 1 (reference)
  #-------------------------------------------------------
  rf1 <- RFsimulate(model, 
                    x = x, 
                    y = y, 
                    T = t,
                    printlevel = 0)
  cmean <- 10
  
  #Non-disturbed field
  vec10 <- scale(rf1$variable1) + cmean
  #Creates bidimensional array
  RF10 <- array(vec10,dim = c(xl*yl,tl))
  #Means
  smeans10 <- sapply(1:10,FUN = function(mth) {
    apply(RF10[,seq(from = mth,
                    by = 10,
                    length.out = 10)],1,sd)
  }
  ) %>% as.vector()
  
  #Disturbed field
  rf1d <- exp(1 + 0.5 * rf1$variable1)
  vec11 <- scale(rf1d) + cmean
  #Creates bidimensional array
  RF11 <- array(vec11,dim = c(xl*yl,tl))
  #Means
  smeans11 <- sapply(1:10,FUN = function(mth) {
    apply(RF11[,seq(from = mth,
                    by = 10,
                    length.out = 10)],1,sd)
  }
  ) %>% as.vector()
  
  #-------------------------------------------------------
  #Creates spatio-temporal field 2 (model)
  #-------------------------------------------------------
  
  rf2 <-  RFsimulate(model, 
                     x = x, 
                     y = y, 
                     T = t,
                     printlevel = 0)
  
  #Defines relations betwenn fields (i counter)
  #-------------------------------------------------------
  for(i in 1:k){
  
    print(i)
  #Correlation
  rho <- Comb[i,1]
  #Standard deviation ratio
  lambda <- Comb[i,2]
  #Bias
  epsilon <- Comb[i,3]
  #Disturbed
  dtb <- Comb[i,4] 

  #-------------------------------------------------------
  #Defines correlation with field 1
  #-------------------------------------------------------
  
  #Defines perturbation
  vec1 <- dtb*vec11 + (1-dtb)*vec10 
  
  n <- length(vec1)
  theta <- acos(rho) #angle
  vec0 <- rf2$variable1       
  X <- cbind(vec1, vec0)
  Xctr <- scale(X, 
                 center = TRUE, 
                 scale = FALSE)
  Id <- diag(n) #Identity matrix
  Q <- qr.Q(qr(Xctr[ , 1, drop = FALSE])) #QR decompositions
  P <- tcrossprod(Q)          # = Q Q' projection over vec1
  x2o <- (Id-P) %*% Xctr[ , 2]
  Xc2 <- cbind(Xctr[ , 1], x2o)                
  Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  #scale to 1
  
  x2 <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]
  vec2 <- x2/sd(x2)*lambda - mean(x2) + epsilon + cmean
  #Store models (series)
  mods.s[i,] <- vec2
  # #Creates bidimensional array
  RF2 <- array(vec2,dim = c(xl*yl,tl))
  #Means------------------------------------------------
  #Field2------------------------------------------------
  smeans2 <- sapply(1:10,FUN = function(mth) {
    apply(RF2[,seq(from = mth,
                   by = 10,
                   length.out = 10)],1,sd)
    }
    ) %>% as.vector()
  #Store models (cycles)
  mods.c[i,] <- smeans2

  }
  
  # save(mods.s,file = "mods.R")
  # save(mods.c,file = "modc.R")
  #load("mods.R")
  #load("modc.R")
 
  #-------------------------------------------------------
  #Chart arguments
  #-------------------------------------------------------  
   
  col.lab <- c(
    `-0.9` = "rho == -0.9",
    `-0.3` = "rho == -0.3",
    `0.3` = "rho == +0.3",
    `0.9` = "rho == +0.9"
  )
  
  row.lab <- c(
    `0` = "delta == 0.0",
    `1` = "delta == 1.0",
    `5` = "delta == 5.0"
  )
  
  #Skill Score
  dist0 <- sapply(1:(k/2),function(i) sscore(vec10,mods.s[i,]))
  dist1 <- sapply((k/2 + 1):k,function(i) sscore(vec11,mods.s[i,]))
  dist <- c(dist0,dist1) %>% as.data.frame()
  dist <- apply(dist,1,deltaf)
  M <- cbind(Comb,dist)
  fctplot(M,"SS")
  ggsave("sscore.png",
         dpi = "retina",
         units = "cm",
         width = 12,
         height = 10)
  
  #Spaef
  dist0 <- sapply(1:(k/2),function(i) spaef(vec10,mods.s[i,]))
  dist1 <- sapply((k/2 + 1):k,function(i) spaef(vec11,mods.s[i,]))
  dist <- c(dist0,dist1) %>% as.data.frame()
  dist <- apply(dist,1,deltaf)
  M <- cbind(Comb,dist)
  fctplot(M,"SPAEF")
  ggsave("spaef.png",
         dpi = "retina",
         units = "cm",
         width = 12,
         height = 10)
  
  #Wasserstein
  dist0 <- sapply(1:(k/2),function(i) wassd(vec10,mods.s[i,]))
  dist1 <- sapply((k/2 + 1):k,function(i) wassd(vec11,mods.s[i,]))
  dist <- c(dist0,dist1) %>% as.data.frame()
  dist <- apply(dist,1,deltaf)
  M <- cbind(Comb,dist)
  fctplot(M,"WD")
  ggsave("wass.png",
         dpi = "retina",
         units = "cm",
         width = 12,
         height = 10)
  
  #Spaef/Wasserstein
  dist0 <- sapply(1:(k/2),function(i) waspaef(vec10,mods.s[i,]))
  dist1 <- sapply((k/2 + 1):k,function(i) waspaef(vec11,mods.s[i,]))
  dist <- c(dist0,dist1) %>% as.data.frame()
  dist <- apply(dist,1,deltaf)
  M <- cbind(Comb,dist)
  fctplot(M,"WASPAEF")
  ggsave("spaefws.png",
         dpi = "retina",
         units = "cm",
         width = 12,
         height = 10)
  
  #Functional analysis
  #Undisturbed field
  #Field1
  #matplot(t(RF1),type = "l")
  #Transform to fdata
  fs10 <- na.omit(fdata(RF10))
  #Calculation of optimal quantity of base function
  nb10 <- optim.basis(fs10,type.basis = "bspline")$numbasis.opt
  #Smooths functions
  fs10.smth.d <- d2f.bspline(df = t(RF10),
                             nb = nb10,
                             av = t,
                             rv = c(1,tl))
  #Field2
  #matplot(t(RF2),type = "l")
  ##Transform to fdata
  dist10 <- sapply(1:(k/2),FUN = function(i){
    RF2 <- array(mods.s[i,],dim = c(xl*yl,tl))
    fs2 <- na.omit(fdata(RF2))
    #Calculation of optimal quantity of base function
    nb2 <- optim.basis(fs2,type.basis = "bspline")$numbasis.opt
    #Smooths functions
    fs2.smth.d <- d2f.bspline(df = t(RF2),
                               nb = nb2,
                               av = t,
                               rv = c(1,tl))
    
    ks1 <- kdpr(t(fs10.smth.d),t(fs2.smth.d),1,I)
    ks2 <- kdsr(fs10.smth.d,fs2.smth.d,"euclidean")
      
    return(list(ks1,ks2))
    
  })
  #Disturbed field
  #Field1
  #matplot(t(RF1),type = "l")
  #Transform to fdata
  fs11 <- na.omit(fdata(RF11))
  #Calculation of optimal quantity of base function
  nb11 <- optim.basis(fs11,type.basis = "bspline")$numbasis.opt
  #Smooths data
  fs11.smth.d <- d2f.bspline(df = t(RF11),
                            nb = nb11,
                            av = t,
                            rv = c(1,tl))
  
  #Field2
  #matplot(t(RF2),type = "l")
  #Transform to fdata
  dist11 <- sapply((k/2 + 1):k,FUN = function(i){
    RF2 <- array(mods.s[i,],dim = c(xl*yl,tl))
    fs2 <- na.omit(fdata(RF2))
    #Calculation of optimal quantity of base function
    nb2 <- optim.basis(fs2,type.basis = "bspline")$numbasis.opt
    #Smooths data
    fs2.smth.d <- d2f.bspline(df = t(RF2),
                              nb = nb2,
                              av = t,
                              rv = c(1,tl))
    
    ks1 <- kdpr(t(fs11.smth.d),t(fs2.smth.d),1,I)
    ks2 <- kdsr(fs11.smth.d,fs2.smth.d,"euclidean")
    
    return(list(ks1,ks2))
  })
  
  #KDPR
  dist <- c(unlist(dist10[1,]),unlist(dist11[1,]))
  M <- cbind(Comb,dist)
  fctplot(M,"KDPR")
  ggsave("ksp.png",
         dpi = "retina",
         units = "cm",
         width = 12,
         height = 10)
  #KDSR
  dist <- c(unlist(dist10[2,]),unlist(dist11[2,]))
  M <- cbind(Comb,dist)
  fctplot(M,"KDSR")
  ggsave("ksd.png",
         dpi = "retina",
         units = "cm",
         width = 12,
         height = 10)
  
 