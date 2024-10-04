#Functions

deltaf <- function(x){
  sqrt(sum((1-x)^2))
}

sscore <- function(o,m){
  c1 <- cor(o,m,use = "pairwise.complete.obs")
  o <- na.omit(o)
  m <- na.omit(m)
  c2 <- c1-(sd(m)/sd(o))
  c3 <- (mean(m)-mean(o))/sd(o)
  ss <- c1^2 - c2^2 - c3^2
  return(ss)
}


spaef <- function(o,m){
  c1 <- 1*(cor(o,m,use = "pairwise.complete.obs")-1)
  o <- na.omit(o)
  m <- na.omit(m)
  c2 <- 1*(((sd(m)/mean(m))/(sd(o)/mean(o)))-1)
  h1 <- hist(scale(o),
             plot = F)
  h2 <- hist(scale(m),
             plot = F)
  nbins <- nclass.FD(scale(m))
  brmax <- max(h1$breaks,h2$breaks)
  brmin <- min(h1$breaks,h2$breaks)
  br <- seq(brmin,brmax,length.out = nbins)
  
  h1 <- hist(scale(o),
             breaks = br,
             plot = F)
  h2 <- hist(scale(m),
             breaks = br,
             plot = F)
  
  c3 <- intersect.dist(h1,h2)
  sf <- 1-sqrt(c1^2+c2^2+c3^2)
  return(sf)
}

wassd <- function(o,m){
  o <- na.omit(o)
  m <- na.omit(m)
  delta <- wasserstein1d(o,m,p = 2)
  return(1-delta)
}

waspaef <- function(o,m){
  c1 <- 1*(cor(o,m,use = "pairwise.complete.obs")-1)
  o <- na.omit(o)
  m <- na.omit(m)
  c2 <- 1*((sd(m)/sd(o))-1)
  so <- scale(o) + mean(o)
  sm <- scale(m) + mean(m)
  c3 <- wasserstein1d(so,sm,p = 2)
  
  wa <- 1-sqrt(c1^2+c2^2+c3^2)

  return(wa)
}

d2f.bspline <- function(df,nb,rv,av){
  basis.obj <- create.bspline.basis(rangeval=rv,
                                      nbasis=nb)
  smooth.obj <- smooth.basis.sparse.mod(argvals = av,
                                    y = df,
                                    fdParobj = basis.obj)
  eval.d <- eval.fd(av, smooth.obj) %>% t()
  return(eval.d)
}

depth <-  function(g, fmat) {
  
  fn <-  ncol(fmat)
  depth <-  rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff <-  abs(sum(sign(g[row] - fmat[row,])))
    depth[row] <-  1 - (diff / fn)
  }
  
  return(depth)
}

int_depth <-  function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}

kdpr <-  function(fmat, gmat,e = 1,fn = I) {
  ff.xd <-  int_depth(fmat, fmat)
  fg.xd <-  int_depth(fmat, gmat)
  
  gg.xd <-  int_depth(gmat, gmat)
  gf.xd <-  int_depth(gmat, fmat)
  
  ff.cdf <-  sapply(ff.xd, function(y) mean(fn(ff.xd - y)^e))
  gf.cdf <-  sapply(ff.xd, function(y) mean(fn(gf.xd - y)^e))
  fg.cdf <-  sapply(gg.xd, function(y) mean(fn(fg.xd - y)^e))
  gg.cdf <-  sapply(gg.xd, function(y) mean(fn(gg.xd - y)^e))
  
  ksf <-  sum(abs(ff.cdf - gf.cdf))
  ksg <-  sum(abs(gg.cdf - fg.cdf))
  kd <- max(ksf,ksg)
  
  return(kd)
}

kdsr <- function(fmat,gmat,ds = "euclidean"){
  
  ff.cdf <- rowMeans(metric.dist(fmat,method = ds))
  gg.cdf <- rowMeans(metric.dist(gmat,method = ds))
  fg.cdf <- rowMeans(metric.dist(fmat,gmat,method = ds))
  
  ksf <-  sum(abs(ff.cdf - fg.cdf))
  ksg <-  sum(abs(gg.cdf - fg.cdf))
  kd <- max(ksf,ksg)

  return(kd)
}


smooth.basis.sparse.mod <- function(argvals, y, fdParobj, fdnames=NULL, covariates=NULL, 
                                    method="qr", dfscale=1 ){
  
  if (is.fdPar(fdParobj)) {
    basisobj = fdParobj$fd$basis
  } else {
    if (is.fd(fdParobj)) {
      basisobj = fdParobj$basis
    } else {
      if (is.basis(fdParobj)) {
        basisobj = fdParobj
      } else {
        stop("fdParobj is not a fdPar, fd, or a basis object.")
      }
    }
  }
  coefs = matrix(0, nrow = basisobj$nbasis, ncol = dim(y)[2])
  for(i in 1:dim(y)[2]){
    curve = y[,i]
    curve.smooth = smooth.basis(argvals[!is.na(curve)],curve[!is.na(curve)],
                                basisobj, covariates, method)
    coefs[,i] = curve.smooth$fd$coefs
  }
  datafd = fd(coefs,basisobj, fdnames)
  
  return(datafd)
}

fctplot <- function(A,mname){
  As <- A %>%
    mutate(avg = rowMeans(dplyr::select(.,5:5))) %>%
    dplyr::select(rho,lambda,epsilon,dbt,avg) %>%
    mutate_at(vars(-avg),factor)
  bc <- ggplot(As, aes(x=lambda, y=avg, group = dbt, color = dbt)) + 
    geom_line(linetype = "dotted")+
    geom_point(aes(shape=dbt),
               size = 1.2) +
    labs(y = mname, x = expression(lambda)) +
    scale_color_manual(labels = c("Undisturbed", "Disturbed"),
                       values = c("#440154","#21918c")) +
    scale_shape_manual(labels = c("Undisturbed", "Disturbed"),
                       values = c(16,17)) +
    facet_grid(epsilon ~ rho,
               scales = "free_y",
               labeller = labeller(rho = as_labeller(col.lab,
                                                     label_parsed),
                                   epsilon = as_labeller(row.lab,
                                                         label_parsed))) +
    theme_bw() +
    theme(axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6),
          strip.text.y.right = element_text(angle = 90),
          legend.title = element_blank(),
          legend.position="bottom",
          legend.margin=margin(t=-10))
  
  return(list(data = As,plot = bc))
}

smooth.pos.mod <- function (argvals, y, WfdParobj, wtvec = rep(1, n), conv = 1e-04, 
                            iterlim = 50, dbglev = 1) 
{
  if (!is.numeric(argvals)) 
    stop("ARGVALS is not numeric.")
  argvals <- as.vector(argvals)
  if (length(argvals) < 2) 
    stop("ARGVALS does not contain at least two values.")
  n <- length(argvals)
  onesobs <- matrix(1, n, 1)
  if (n < 2) 
    stop("ARGVALS does not contain at least two values.")
  y = as.matrix(y)
  ychk <- ycheck(y, n)
  y <- ychk$y
  ncurve <- ychk$ncurve
  nvar <- ychk$nvar
  ndim <- ychk$ndim
  WfdParobj <- fdParcheck(WfdParobj, curve)
  lambda <- WfdParobj$lambda
  Wfdobj <- WfdParobj$fd
  Lfdobj <- WfdParobj$Lfd
  basisobj <- Wfdobj$basis
  nbasis <- basisobj$nbasis
  wtvec <- wtcheck(n, wtvec)$wtvec
  climit <- c(rep(-400, nbasis), rep(400, nbasis))
  coef0 <- matrix(0,nbasis,ncurve)
  active <- 1:nbasis
  if (lambda > 0) 
    Kmat <- lambda * eval.penalty(basisobj, Lfdobj)
  else Kmat <- matrix(0, nbasis, nbasis)
  if (ndim == 2) {
    coef <- matrix(0, nbasis, ncurve)
  }
  else {
    coef <- array(0, c(nbasis, ncurve, nvar))
  }
  if (ncurve > 1 || nvar > 1) 
    Flist <- vector("list", ncurve * nvar)
  else Flist <- NULL
  for (ivar in 1:nvar) {
    if (ndim == 2) {
      sclfac <- mean(c(y)^2)
    }
    else {
      sclfac <- mean(c(y[, , ivar])^2)
    }
    for (icurve in 1:ncurve) {
      if (ndim == 2) {
        yi <- y[, icurve]
        cveci <- coef0[, icurve]
      }
      else {
        yi <- y[, icurve, ivar]
        cveci <- coef0[, icurve, ivar]
      }
      result <- PENSSEfun(argvals, yi, basisobj, cveci, 
                          Kmat, wtvec)
      PENSSE <- result[[1]]
      DPENSSE <- result[[2]]
      f0 <- PENSSE
      gvec0 <- DPENSSE
      Foldlist <- list(f = f0, grad = gvec0, norm = sqrt(mean(gvec0^2)))
      hmat0 <- PENSSEhess(argvals, yi, basisobj, cveci, 
                          Kmat, wtvec)
      deltac <- -solve(hmat0, gvec0)
      cosangle <- -sum(gvec0 * deltac)/sqrt(sum(gvec0^2) * 
                                              sum(deltac^2))
      iternum <- 0
      status <- c(iternum, Foldlist$f, -PENSSE, Foldlist$norm)
      if (dbglev >= 1) {
        if (ncurve > 1 || nvar > 1) {
          if (ncurve > 1 && nvar > 1) {
            cat("\n")
            curvetitle <- paste("Results for curve", 
                                icurve, "and variable", ivar)
          }
          if (ncurve > 1 && nvar == 1) {
            cat("\n")
            curvetitle <- paste("Results for curve", 
                                icurve)
          }
          if (ncurve == 1 && nvar > 1) {
            cat("\n")
            curvetitle <- paste("Results for variable", 
                                ivar)
          }
          cat("\n")
          cat(curvetitle)
        }
        cat("\n")
        cat("\nIter.   PENSSE   Grad Length")
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2], 4))
        cat("      ")
        cat(round(status[3], 4))
      }
      MAXSTEPITER <- 10
      MAXSTEP <- 400
      trial <- 1
      linemat <- matrix(0, 3, 5)
      Flisti <- Foldlist
      gvec <- gvec0
      dbgwrd <- dbglev > 1
      if (iterlim == 0) {
        cat("\n")
      }
      else {
        for (iter in 1:iterlim) {
          iternum <- iternum + 1
          dblwrd <- rep(FALSE, 2)
          limwrd <- rep(FALSE, 2)
          stpwrd <- FALSE
          ind <- 0
          ips <- 0
          linemat[2, 1] <- sum(deltac * Flisti$grad)
          sdg <- sqrt(sum(deltac^2))
          deltac <- deltac/sdg
          dgsum <- sum(deltac)
          linemat[2, 1] <- linemat[2, 1]/sdg
          linemat[, 1:4] <- outer(c(0, linemat[2, 1], 
                                    Flisti$f), rep(1, 4))
          stepiter <- 0
          if (dbglev >= 2) {
            cat("\n")
            cat(paste("                 ", stepiter, 
                      "  "))
            cat(format(round(t(linemat[, 1]), 6)))
          }
          if (linemat[2, 1] >= 0) {
            print("Initial slope nonnegative.")
            ind <- 3
            break
          }
          if (linemat[2, 1] >= -1e-05) {
            if (dbglev > 1) 
              print("Initial slope too small")
            break
          }
          linemat[1, 5] <- trial
          for (stepiter in 1:MAXSTEPITER) {
            limflg <- 0
            result <- stepchk(linemat[1, 5], cveci, deltac, 
                              limwrd, ind, climit, active, dbgwrd)
            linemat[1, 5] <- result[[1]]
            ind <- result[[2]]
            limwrd <- result[[3]]
            if (linemat[1, 5] <= 1e-07) {
              Flisti <- Foldlist
              cvecnew <- cveci
              gvecnew <- gvec
              if (dbglev >= 2) {
                print("Stepsize too small")
                print(linemat[1, 5])
              }
              if (limflg) 
                ind <- 1
              else ind <- 4
              break
            }
            cvecnew <- cveci + linemat[1, 5] * deltac
            result <- PENSSEfun(argvals, yi, basisobj, 
                                cvecnew, Kmat, wtvec)
            PENSSE <- result[[1]]
            DPENSSE <- result[[2]]
            Flisti$f <- PENSSE
            gvecnew <- DPENSSE
            Flisti$grad <- gvecnew
            Flisti$norm <- sqrt(mean(gvecnew^2))
            linemat[3, 5] <- Flisti$f
            linemat[2, 5] <- sum(deltac * gvecnew)
            if (dbglev >= 2) {
              cat("\n")
              cat(paste("                 ", stepiter, 
                        "  "))
              cat(format(round(t(linemat[, 5]), 6)))
            }
            result <- stepit(linemat, ips, dblwrd, MAXSTEP)
            linemat <- result[[1]]
            ips <- result[[2]]
            ind <- result[[3]]
            dblwrd <- result[[4]]
            trial <- linemat[1, 5]
            if (ind == 0 | ind == 5) 
              break
          }
          cveci <- cvecnew
          gvec <- gvecnew
          if (abs(Flisti$f - Foldlist$f) < sclfac * conv) {
            cat("\n")
            break
          }
          if (Flisti$f >= Foldlist$f) 
            break
          hmat <- PENSSEhess(argvals, yi, basisobj, cveci, 
                             Kmat, wtvec)
          deltac <- -solve(hmat, gvec)
          cosangle <- -sum(gvec * deltac)/sqrt(sum(gvec^2) * 
                                                 sum(deltac^2))
          if (cosangle < 0) {
            if (dbglev > 1) 
              print("cos(angle) negative")
            deltac <- -gvec
          }
          Foldlist <- Flisti
          status <- c(iternum, Flisti$f, Flisti$norm)
          if (dbglev >= 1) {
            cat("\n")
            cat(iternum)
            cat("        ")
            cat(round(status[2], 4))
            cat("      ")
            cat(round(status[3], 4))
          }
        }
      }
      if (ndim == 2) {
        coef[, icurve] <- cveci
      }
      else {
        coef[, icurve, ivar] <- cveci
      }
      if (ncurve == 1 && nvar == 1) {
        Flist <- Flisti
      }
      else {
        Flist[[(ivar - 1) * ncurve + icurve]] <- Flisti
      }
    }
  }
  Wfdobj <- fd(coef, basisobj)
  posFd <- list(Wfdobj = Wfdobj, Flist = Flist, argvals = argvals, 
                y = y)
  class(posFd) <- "posfd"
  return(posFd)
}


#  ---------------------------------------------------------------

PENSSEfun <- function(argvals, yi, basisobj, cveci, Kmat, wtvec) {
  #  Computes the log likelihood and its derivative with
  #    respect to the coefficients in CVEC
  n       <- length(argvals)
  nbasis  <- basisobj$nbasis
  phimat  <- getbasismatrix(argvals, basisobj, 0)
  Wvec    <- phimat %*% cveci
  EWvec   <- exp(Wvec)
  res     <- yi - EWvec
  PENSSE  <- mean(wtvec*res^2) + t(cveci) %*% Kmat %*% cveci
  DPENSSE <- -2*crossprod(phimat,wtvec*res*EWvec)/n + 2*Kmat %*% cveci
  return( list(PENSSE, DPENSSE) )
}

#  ---------------------------------------------------------------

PENSSEhess <- function(argvals, yi, basisobj, cveci, Kmat, wtvec) {
  #  Computes the expected Hessian
  n        <- length(argvals)
  nbasis   <- basisobj$nbasis
  phimat   <- getbasismatrix(argvals, basisobj, 0)
  Wvec     <- phimat %*% cveci
  EWvec    <- exp(Wvec)
  D2PENSSE <- 2*t(phimat) %*% diag(as.numeric(wtvec*EWvec^2)) %*% phimat/n + 2*Kmat
  return(D2PENSSE)
}

process.raster <- function(f,v,d1,d2,lmts,lp,tmpl){
  b <- brick(f, varname = v, lvar = 3, level = 1) #Creates raster object
  tid <- which(getZ(b) >= d1 & getZ(b) <  d2)  #Generates time index
  btsr <- b %>% 
    raster::subset(tid) %>% #Time filter
    raster::crop(lmts) %>% #Cuts area of interest
    '*'(lp[1]) %>%
    '+'(lp[2]) %>%
    raster::resample(tmpl,method="bilinear")  #resample
  btsr@z$Date <- b@z$Date[tid]
  return(btsr)
}

group.raster <- function(btsr,ind,aggfun){
  aggb <- stackApply(btsr,
                     indices = ind, 
                     fun = aggfun,
                     na.rm = TRUE)
  return(aggb)
}

# Generates teleconnection patterns
tlc.ptn <- function(ts,tsm,sls,yrq){
  qrt.seq <- rep(1:4,yrq)
  tlc.matrix <- matrix(NA,
                       nrow = sls,
                       ncol = 4)
  for(qrt in 1:4){
    station.ind <- which(qrt.seq==qrt)
    for(i in 1:sls){
      tlc.matrix[i,qrt] <- cor(ts[station.ind],
                               tsm[i,station.ind])
    }
  }
  return(tlc.matrix)
}

#Calculates different measures over the list
calc.measure <- function(olist,mlist,fn){
  colq <- ncol(olist)
  results <- matrix(NA,nrow = mdq,ncol = colq)
  
  for(i in 1:mdq){
    v <- sapply(1:colq,FUN = function(j){
      o <- olist[,j]
      m <- mlist[[i]][,j]
      fn(o,m)
    }
    )
    results[i,] <- v
  }
  return(results)
}

#Euclidean Norm
normvec <- function(x) sqrt(sum(x^2))

#Boxplot per month or season
bxp <- function(df,clnames,hglmds,lyt){
  cycledf <- data.frame(df) %>%
    `colnames<-`(clnames) %>% 
    bind_cols(namesdf) %>% 
    mutate(nmod = factor(1:32)) %>%
    #arrange(name,variant) %>%
    mutate(ordmod = factor(1:32)) %>% 
    mutate(identifier = str_c(ordmod,".",name,"\n",variant)) %>%
    mutate(identifier = factor(identifier,
                               levels = identifier)) %>%
    pivot_longer(!c(name,variant,nmod,ordmod,identifier),
                 names_to = "timeu",
                 values_to = "waspaef") %>%
    mutate(timeu = factor(timeu,
                          levels = clnames))
  
  highmod <- cycledf %>%
    filter(nmod %in% hglmds)
  
  gg <- ggplot(cycledf, 
               aes(x=timeu, 
                   y=waspaef)) + 
    geom_boxplot() +
    geom_point(data = highmod,
               aes(x=timeu, 
                   y=waspaef,
                   color = identifier,
                   shape = identifier)) +
    scale_color_viridis(discrete=TRUE) +
    scale_shape_manual(values=c(3,4,7,8,9,10)) +
    labs(x = NULL,
         y = "WASPAEF",
         color = NULL,
         shape = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.text=element_text(size=7),
          plot.margin = margin(2, 2, 2, 2)) +
    guides(color = guide_legend(nrow = lyt[1],ncol = lyt[2], byrow = T),
           shape = guide_legend(nrow = lyt[1],ncol = lyt[2], byrow = T))
  return(gg)
}

