#Analytic Code for "Healthcare Costs of Treating Privately Insured Patients with Acute Myeloid Leukemia in the United States from 2004 to 2014: A Generalized Additive Modeling Approach"
#gam reference: https://cran.r-project.org/web/packages/gam/gam.pdf 

install.packages("gam")
rm(list=ls())
library("gam")
versnum <- "ALL_20190311-1" #Analyzes HCT and chemo only patients

ffl <- function(vec, first, last) {
  aa1 <- match(vec, table = last)
  aa1[is.na(aa1)] <- 0
  aa2 <- cumsum(aa1)
  l <- which.max(match(aa2, 2))
  f <- l + first
  return(c(f,l))
}

hct.store.eq <- list()
chemo.store.eq <- list()
hct.store.plots <- list()
chemo.store.plots <- list()

#set working directory
setwd ("C:/gam")

#read in data : Note: anaytic dataset not provided due to contractual limitations
chemo<-read.csv("Chemocumulativecostsanalytic.csv")
hct<-read.csv("HCTcumulativecostsanalytic.csv")
chemo.df <- data.frame(chemo$Patid, chemo$X1_chemo, chemo$X2_HCT, chemo$X_death, chemo$X_enroll, chemo$deathleLDE, chemo$cum_sum)
names(chemo.df) <- c("Patid","X1_chemo","X2_HCT","X_death", "X_enroll", "deathleLDE", "cum_sum")
hct.df <- data.frame(hct$Patid,hct$X1_chemo, hct$X2_HCT, hct$X_death, hct$X_enroll, hct$deathleLDE, hct$cum_sum)
names(hct.df) <- c("Patid","X1_chemo","X2_HCT","X_death", "X_enroll", "deathleLDE", "cum_sum")

#HCT 

hct.deathna <- hct.df$X_death
hct.deathna[is.na(hct.deathna)]<- -30
hct.deathna[hct.deathna < -30] <- -30
hct.tHCT <- hct.df$X2_HCT
hct.tHCT[hct.tHCT < -30] <- -30
hct.tLDE <- hct.df$X_enroll
hct.tLDE[hct.tLDE < -30] <- -30
hct.tLDE[hct.tLDE > 0] <- -30
hct.tLDE[hct.df$deathleLDE==1] <- -30
hct.df$t_death <- hct.deathna
hct.df$t_LDE <- hct.tLDE
hct.df$t_HCT <- hct.tHCT

hct.eq <- gam(cum_sum ~ ns(x = X1_chemo,knots = c(30,60,90,120,150,180,365,550,730,1095),Boundary.knots = c(0,max(X1_chemo)))
              + ns(x = t_HCT,knots = c(-14,0,14,30,100,180,365,550),Boundary.knots = c(-30,2*365))
              + ns(x = t_death,knots = c(-14,-7),Boundary.knots = c(-30,0))
              + ns(x = t_LDE,knots = c(-14,-7),Boundary.knots = c(-30,0)), family=Gamma(link=log), data=hct.df, model = FALSE)

#Chemo Only

chemo.deathna <- chemo.df$X_death
chemo.deathna[is.na(chemo.deathna)]<- -30
chemo.deathna[chemo.deathna < -30] <- -30
chemo.tLDE <- chemo.df$X_enroll
chemo.tLDE[chemo.tLDE < -30] <- -30
chemo.tLDE[chemo.tLDE > 0] <- -30
chemo.tLDE[chemo.df$deathleLDE==1] <- -30
chemo.df$t_death <- chemo.deathna
chemo.df$t_LDE <- chemo.tLDE

chemo.eq <- gam(cum_sum ~ ns(x = X1_chemo,knots = c(30,60,90,120,150,180,365,550,730,1095),Boundary.knots = c(0,max(X1_chemo)))
                + ns(x = t_death,knots = c(-14,-7),Boundary.knots = c(-30,0))
                + ns(x = t_LDE,knots = c(-14,-7),Boundary.knots = c(-30,0)), family=Gamma(link=log), data=chemo.df, model = FALSE)

#Plot GAM results
pdf(file=paste("basecase_results_",versnum,".pdf",sep = ""))
plot(hct.eq,residuals = FALSE,se = TRUE,rugplot = FALSE, scale = 5, ylim = c(-3,2), ylab ="Contribution to log-Cost", main = "HCT" )
plot(chemo.eq,residuals = FALSE,se = TRUE,rugplot = FALSE, scale = 5, ylim = c(-3,2), ylab ="Contribution to log-Cost", main = "Chemo Only" )
dev.off()

hct.bs_plot <- preplot(hct.eq,residuals = FALSE,se = FALSE,rugplot = FALSE, scale = 5, ylim = c(-3,2), ylab ="Contribution to log-Cost")
hct.chemo_x <- hct.bs_plot[[1]]$x
hct.chemo_y <- hct.bs_plot[[1]]$y
hct.hct_x <- hct.bs_plot[[2]]$x
hct.hct_y <- hct.bs_plot[[2]]$y
hct.death_x <- hct.bs_plot[[3]]$x
hct.death_y <- hct.bs_plot[[3]]$y
hct.lde_x <- hct.bs_plot[[4]]$x
hct.lde_y <- hct.bs_plot[[4]]$y
hct.store.plots[[1]] <- list(chemo_xy = cbind(hct.chemo_x[(which.max(hct.chemo_x)-max(hct.chemo_x)):which.max(hct.chemo_x)],hct.chemo_y[(which.max(hct.chemo_x)-max(hct.chemo_x)):which.max(hct.chemo_x)]),
                             hct_xy = cbind(hct.hct_x[(which.max(hct.hct_x)-(max(hct.hct_x)-min(hct.hct_x))):which.max(hct.hct_x)], hct.hct_y[(which.max(hct.hct_x)-(max(hct.hct_x)-min(hct.hct_x))):which.max(hct.hct_x)]),
                             death_xy = cbind(hct.death_x[(which.max(hct.death_x)-(max(hct.death_x)-min(hct.death_x))):which.max(hct.death_x)], hct.death_y[(which.max(hct.death_x)-(max(hct.death_x)-min(hct.death_x))):which.max(hct.death_x)]),
                             lde_xy = cbind(hct.lde_x[(which.max(hct.lde_x)-(max(hct.lde_x)-min(hct.lde_x))):which.max(hct.lde_x)], hct.lde_y[(which.max(hct.lde_x)-(max(hct.lde_x)-min(hct.lde_x))):which.max(hct.lde_x)]))

chemo.bs_plot <- preplot(chemo.eq,residuals = FALSE,se = FALSE,rugplot = FALSE, scale = 5, ylim = c(-3,2), ylab ="Contribution to log-Cost")
chemo.chemo_x <- chemo.bs_plot[[1]]$x
chemo.chemo_y <- chemo.bs_plot[[1]]$y
chemo.death_x <- chemo.bs_plot[[2]]$x
chemo.death_y <- chemo.bs_plot[[2]]$y
chemo.lde_x <- chemo.bs_plot[[3]]$x
chemo.lde_y <- chemo.bs_plot[[3]]$y
chemo.store.plots[[1]] <- list(chemo_xy = cbind(chemo.chemo_x[(which.max(chemo.chemo_x)-max(chemo.chemo_x)):which.max(chemo.chemo_x)],chemo.chemo_y[(which.max(chemo.chemo_x)-max(chemo.chemo_x)):which.max(chemo.chemo_x)]),
                               death_xy = cbind(chemo.death_x[(which.max(chemo.death_x)-(max(chemo.death_x)-min(chemo.death_x))):which.max(chemo.death_x)], chemo.death_y[(which.max(chemo.death_x)-(max(chemo.death_x)-min(chemo.death_x))):which.max(chemo.death_x)]),
                               lde_xy = cbind(chemo.lde_x[(which.max(chemo.lde_x)-(max(chemo.lde_x)-min(chemo.lde_x))):which.max(chemo.lde_x)], chemo.lde_y[(which.max(chemo.lde_x)-(max(chemo.lde_x)-min(chemo.lde_x))):which.max(chemo.lde_x)]))

#Generate weekly grid for storing predictions (30 days to 5 years) of cumulative cost as of death date for patients known to have died
#Note: HCT dates after death or LDE are nonsense and will be set to NA ex post
HCTgridval <- seq(30,5*365,7)
deathgridval <- seq(30,5*365,7)
LDEgridval <- seq(30,5*365,7)
x1 <- HCTgridval
x2 <- deathgridval
n1 <- length(x1)
n2 <- length(x2)
x12_t <- matrix(data = NA,nrow = n1*n2,ncol = 2)
for (i in 1:n1){
  x12_t[(n2*(i-1)+1):(n2*(i-1)+n2),] <- cbind(rep(x1[i],n2),x2)
}

#Clear values to save storage space : not necessary for lookup tool
hct.eq$data <- NULL
hct.eq$residuals <- NULL
hct.eq$fitted.values <- NULL
hct.eq$effects <- NULL
hct.eq$weights <- NULL
hct.eq$prior.weights <- NULL
hct.eq$y <- NULL
hct.eq$qr$qr <- NULL
hct.eq$additive.predictors <- NULL
chemo.eq$data <- NULL
chemo.eq$residuals <- NULL
chemo.eq$fitted.values <- NULL
chemo.eq$effects <- NULL
chemo.eq$weights <- NULL
chemo.eq$prior.weights <- NULL
chemo.eq$y <- NULL
chemo.eq$qr$qr <- NULL
chemo.eq$additive.predictors <- NULL

hct.store.eq[[1]] <- hct.eq
chemo.store.eq[[1]] <- chemo.eq

#Bootstrap hct : must always be run after single run equation above because times are defined above

M <- 1000 #Number of bootstrap iterations
hct.bs_ids <- unique(hct.df$Patid)
chemo.bs_ids <- unique(chemo.df$Patid)
hct.n_all <- length(hct.bs_ids)
chemo.n_all <- length(chemo.bs_ids)

set.seed(314159)

tick <- as.numeric(Sys.time())
for (j in 1:M){
  hct.bs_sampids <- sample(x = hct.bs_ids,size = hct.n_all, replace = TRUE)
  chemo.bs_sampids <- sample(x = chemo.bs_ids,size = chemo.n_all, replace = TRUE)
  hct.bs_all <- data.frame()
  chemo.bs_all <- data.frame()
  for (i in 1: hct.n_all){
    hct.bs_all_i <- hct.df[hct.df$Patid==hct.bs_sampids[i],c("Patid","X1_chemo","t_HCT","t_death","t_LDE","cum_sum")]
    hct.bs_all <- rbind(hct.bs_all,hct.bs_all_i)
  }
  
  for (i in 1: chemo.n_all){
    chemo.bs_all_i <- chemo.df[chemo.df$Patid==chemo.bs_sampids[i],c("Patid","X1_chemo","t_death","t_LDE","cum_sum")]
    chemo.bs_all <- rbind(chemo.bs_all,chemo.bs_all_i)
  }
  
  hct.bs_eq <- gam(cum_sum ~ ns(x = X1_chemo,knots = c(30,60,90,120,150,180,365,550,730,1095),Boundary.knots = c(0,max(X1_chemo)))
                   + ns(x = t_HCT,knots = c(-14,0,14,30,100,180,365,550),Boundary.knots = c(-30,730))
                   + ns(x = t_death,knots = c(-14,-7),Boundary.knots = c(-30,0))
                   + ns(x = t_LDE,knots = c(-14,-7),Boundary.knots = c(-30,0)), family=Gamma(link=log), data=hct.bs_all, model = FALSE)
  
  chemo.bs_eq <- gam(cum_sum ~ ns(x = X1_chemo,knots = c(30,60,90,120,150,180,365,550,730,1095),Boundary.knots = c(0,max(X1_chemo)))
                     + ns(x = t_death,knots = c(-14,-7),Boundary.knots = c(-30,0))
                     + ns(x = t_LDE,knots = c(-14,-7),Boundary.knots = c(-30,0)), family=Gamma(link=log), data=chemo.bs_all, model = FALSE)
  
  hct.bs_plot <- preplot(hct.bs_eq,residuals = FALSE,se = FALSE,rugplot = FALSE, scale = 5, ylim = c(-3,2), ylab ="Contribution to log-Cost")
  hct.chemo_x <- hct.bs_plot[[1]]$x
  hct.chemo_y <- hct.bs_plot[[1]]$y
  hct.hct_x <- hct.bs_plot[[2]]$x
  hct.hct_y <- hct.bs_plot[[2]]$y
  hct.death_x <- hct.bs_plot[[3]]$x
  hct.death_y <- hct.bs_plot[[3]]$y
  hct.lde_x <- hct.bs_plot[[4]]$x
  hct.lde_y <- hct.bs_plot[[4]]$y
  hct.first.chemo.x <- (which.max(hct.chemo_x)-max(hct.chemo_x))
  hct.last.chemo.x <- which.max(hct.chemo_x)
  hct.first.hct.x <- (which.max(hct.hct_x)-max(hct.hct_x))
  hct.last.hct.x <- which.max(hct.hct_x)
  fl_5 <- ffl(hct.death_x,-30,0)
  hct.first.death.x <- fl_5[1]
  hct.last.death.x <- fl_5[2]
  fl_6 <- ffl(hct.lde_x,-30,0)
  hct.first.lde.x <- fl_6[1]
  hct.last.lde.x <- fl_6[2]
  
  hct.store.plots[[j+1]] <- list(chemo_xy = cbind(hct.chemo_x[hct.first.chemo.x:hct.last.chemo.x],hct.chemo_y[hct.first.chemo.x:hct.last.chemo.x]),
                                 hct_xy = cbind(hct.hct_x[hct.first.hct.x:hct.last.hct.x], hct.hct_y[hct.first.hct.x:hct.last.hct.x]),
                                 death_xy = cbind(hct.death_x[hct.first.death.x:hct.last.death.x], hct.death_y[hct.first.death.x:hct.last.death.x]),
                                 lde_xy = cbind(hct.lde_x[hct.first.lde.x:hct.last.lde.x], hct.lde_y[hct.first.lde.x:hct.last.lde.x]))
  
  chemo.bs_plot <- preplot(chemo.bs_eq,residuals = FALSE,se = FALSE,rugplot = FALSE, scale = 5, ylim = c(-3,2), ylab ="Contribution to log-Cost")
  chemo.chemo_x <- chemo.bs_plot[[1]]$x
  chemo.chemo_y <- chemo.bs_plot[[1]]$y
  chemo.death_x <- chemo.bs_plot[[2]]$x
  chemo.death_y <- chemo.bs_plot[[2]]$y
  chemo.lde_x <- chemo.bs_plot[[3]]$x
  chemo.lde_y <- chemo.bs_plot[[3]]$y
  chemo.first.chemo.x <- (which.max(chemo.chemo_x)-max(chemo.chemo_x))
  chemo.last.chemo.x <- which.max(chemo.chemo_x)
  fl_1 <- ffl(chemo.death_x,-30,0)
  chemo.first.death.x <- fl_1[1]
  chemo.last.death.x <- fl_1[2]
  fl_2 <- ffl(chemo.lde_x,-30,0)
  chemo.first.lde.x <- fl_2[1]
  chemo.last.lde.x <- fl_2[2]
  
  chemo.store.plots[[j+1]] <- list(chemo_xy = cbind(chemo.chemo_x[chemo.first.chemo.x:chemo.last.chemo.x],chemo.chemo_y[chemo.first.chemo.x:chemo.last.chemo.x]),
                                   death_xy = cbind(chemo.death_x[chemo.first.death.x:chemo.last.death.x], chemo.death_y[chemo.first.death.x:chemo.last.death.x]),
                                   lde_xy = cbind(chemo.lde_x[chemo.first.lde.x:chemo.last.lde.x], chemo.lde_y[chemo.first.lde.x:chemo.last.lde.x]))
  
  hct.bs_eq$data <- NULL
  hct.bs_eq$residuals <- NULL
  hct.bs_eq$fitted.values <- NULL
  hct.bs_eq$effects <- NULL
  hct.bs_eq$weights <- NULL
  hct.bs_eq$prior.weights <- NULL
  hct.bs_eq$y <- NULL
  hct.bs_eq$qr$qr <- NULL
  hct.bs_eq$additive.predictors <- NULL
  chemo.bs_eq$data <- NULL
  chemo.bs_eq$residuals <- NULL
  chemo.bs_eq$fitted.values <- NULL
  chemo.bs_eq$effects <- NULL
  chemo.bs_eq$weights <- NULL
  chemo.bs_eq$prior.weights <- NULL
  chemo.bs_eq$y <- NULL
  chemo.bs_eq$qr$qr <- NULL
  chemo.bs_eq$additive.predictors <- NULL
  
  hct.store.eq[[j+1]] <- hct.bs_eq
  chemo.store.eq[[j+1]] <- chemo.bs_eq
  
  tock <- as.numeric(Sys.time())
  elapsed.time <- tock - tick
  time.left <- (elapsed.time / j)*(M-j)
  print(paste("Estimated total run time left:", time.left/60, " minutes."))
  flush.console()
}

save(hct.store.eq,hct.store.plots,chemo.store.eq,chemo.store.plotsfile = paste("results_",versnum,".Rdata",sep = ""))

pdf(file=paste("hct.bs_plots_",versnum,".pdf",sep = ""),width=10,height=7.5)
for (k in 1:(M+1)){
  if (k==1){
    plot(hct.store.plots[[k]]$chemo_xy[,1],hct.store.plots[[k]]$chemo_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(hct.store.plots[[k]]$chemo_xy[,1],hct.store.plots[[k]]$chemo_xy[,2], type = "l", col = k)
  }
}
for (k in 1:(M+1)){
  if (k==1){
    plot(hct.store.plots[[k]]$hct_xy[,1],hct.store.plots[[k]]$hct_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(hct.store.plots[[k]]$hct_xy[,1],hct.store.plots[[k]]$hct_xy[,2], type = "l", col = k)
  }
}
for (k in 1:(M+1)){
  if (k==1){
    plot(hct.store.plots[[k]]$death_xy[,1],hct.store.plots[[k]]$death_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(hct.store.plots[[k]]$death_xy[,1],hct.store.plots[[k]]$death_xy[,2], type = "l", col = k)
  }
}
for (k in 1:(M+1)){
  if (k==1){
    plot(hct.store.plots[[k]]$lde_xy[,1],hct.store.plots[[k]]$lde_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(hct.store.plots[[k]]$lde_xy[,1],hct.store.plots[[k]]$lde_xy[,2], type = "l", col = k)
  }
}
dev.off()

pdf(file=paste("chemo.bs_plots_",versnum,".pdf",sep = ""),width=10,height=7.5)
for (k in 1:(M+1)){
  if (k==1){
    plot(chemo.store.plots[[k]]$chemo_xy[,1],chemo.store.plots[[k]]$chemo_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(chemo.store.plots[[k]]$chemo_xy[,1],chemo.store.plots[[k]]$chemo_xy[,2], type = "l", col = k)
  }
}
for (k in 1:(M+1)){
  if (k==1){
    plot(chemo.store.plots[[k]]$death_xy[,1],chemo.store.plots[[k]]$death_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(chemo.store.plots[[k]]$death_xy[,1],chemo.store.plots[[k]]$death_xy[,2], type = "l", col = k)
  }
}
for (k in 1:(M+1)){
  if (k==1){
    plot(chemo.store.plots[[k]]$lde_xy[,1],chemo.store.plots[[k]]$lde_xy[,2], type = "l", ylim = c(-3,2), lwd = 3, col = "black")
  } else {
    lines(chemo.store.plots[[k]]$lde_xy[,1],chemo.store.plots[[k]]$lde_xy[,2], type = "l", col = k)
  }
}
dev.off()

gam_lookup <- function(DATASET,DATE_HCT, DATE_DEATH, DATE_LDE,PLOT,yoffset,xoffset){
  if (DATASET=="HCT"){
    store.eq <- "hct.store.eq"
    maintext <- "HCT Recipient Cohort Expected Cost\n"
  } else {
    store.eq <- "all.store.eq"
    maintext <- "Propensity-Score Weighted Combined Cohort Expected Cost\n"    
  } 
  if(is.na(DATE_LDE)){
    #Death date known
    if (DATE_DEATH >= 365){
      if(DATE_DEATH%%365==0) {
        chemoseq <- seq(from = 365,to = 365*floor(DATE_DEATH/365),by = 365)
      } else {chemoseq <- c(seq(from = 365,to = 365*floor(DATE_DEATH/365),by = 365),DATE_DEATH)}
    } else {chemoseq <- DATE_DEATH}
    deathseq <- chemoseq - DATE_DEATH
    deathseq[deathseq < -30] <- -30
    ldeseq <- rep(-30,length(chemoseq))
  } else {
    if (DATE_LDE >= 365){
      if(DATE_LDE%%365==0) {
        chemoseq <- seq(from = 365,to = 365*floor(DATE_LDE/365),by = 365)
      } else {chemoseq <- c(seq(from = 365,to = 365*floor(DATE_LDE/365),by = 365),DATE_LDE)}
    } else {chemoseq <- DATE_LDE}
    ldeseq <- chemoseq - DATE_LDE
    ldeseq[ldeseq < -30] <- -30
    deathseq <- rep(-30,length(chemoseq))
  }
  if (is.na(DATE_HCT)) {
    hctseq <- rep(-30,length(chemoseq))
  } else{
    hctseq <- chemoseq - DATE_HCT
    hctseq[hctseq < -30] <- -30
  }
  working.dat <- data.frame(X1_chemo = chemoseq, t_HCT = hctseq, t_death = deathseq, t_LDE = ldeseq)
  pred.res<-lapply(get(store.eq), FUN=predict, newdata=working.dat, type = "response")
  hct.bs0 <- do.call(rbind,pred.res[2:1001])
  hct.bs <- matrix(NA,nrow = 1000,ncol=dim(hct.bs0)[2])
  hct.bs[,1] <- hct.bs0[,1]
  for (k in 2:dim(hct.bs0)[2]){
    hct.bs[,k] <- pmax(hct.bs0[,k],hct.bs[,k-1])
  }
  hct.bs.diff <- cbind(hct.bs[,1],t(diff(t(hct.bs))))
  hct.bc <- unlist(pred.res[1])
  hct.bc.diff <- c(hct.bc[1],diff(hct.bc))
  hct.025 <- apply(hct.bs,MARGIN = 2, FUN=quantile, probs = .025)
  hct.975 <- apply(hct.bs,MARGIN = 2, FUN=quantile, probs = .975)
  hct.diff.025 <- apply(hct.bs.diff,MARGIN = 2, FUN=quantile, probs = .025)
  hct.diff.975 <- apply(hct.bs.diff,MARGIN = 2, FUN=quantile, probs = .975)
  if(is.na(DATE_HCT)){
    if (is.na(DATE_DEATH)){
      subtext <- paste("No HCT, Lost to Follow-up on Day", DATE_LDE) 
    } else {
      subtext <- paste("No HCT, Died on Day", DATE_DEATH)
    } 
  } else {
    if (is.na(DATE_DEATH)){
      subtext <- paste("HCT on Day", DATE_HCT, "/ Lost to Follow-up on Day", DATE_LDE)
    } else {
      subtext <- paste("HCT on Day", DATE_HCT, "/ Died on Day", DATE_DEATH)
    }
  }
  if (PLOT==TRUE){
    plot(c(0,chemoseq),c(0,hct.bc),pch=19,axes = FALSE, ylim=c(0,2000000),xlim=c(0,5*365+200),xlab = "Days since Chemotherapy Initiation",ylab = "Cumulative Cost",main = paste(maintext,subtext,sep=""))
    arrows(c(chemoseq), c(hct.025), c(chemoseq), c(hct.975), length=0.05, angle=90, code=3)
    text(x = chemoseq+90+xoffset,y = hct.bc+yoffset,labels = paste("$",format(round(hct.bc,-3),scientific = FALSE,big.mark = ",")),cex=.7, pos=1, offset=.5)
    axis(side = 1, at = seq(0,365*5,365))
    axis(side=2,at=seq(from=0,to=2000000,by=200000), labels = c(format(seq(from=0,to=2000000,by=200000),scientific = FALSE)))
  } 
  return(data.frame(chemoseq,hct.bc,hct.025,hct.975,hct.bc.diff,hct.diff.025,hct.diff.975))
}

# HCT
manuscript_HCT = "C:/  HCT.pdf"
pdf(file= HCT)
par(mfrow=c(3,2), cex=.5)
gam_lookup(DATASET = "HCT",DATE_HCT = 60,DATE_DEATH = 1825,DATE_LDE = NA,PLOT=TRUE,yoffset= -60000,xoffset=50)
gam_lookup(DATASET = "HCT",DATE_HCT = 60,DATE_DEATH = NA,DATE_LDE = 1825,PLOT=TRUE,yoffset= -60000,xoffset = 50)
gam_lookup(DATASET = "HCT",DATE_HCT = 425,DATE_DEATH = 1825,DATE_LDE = NA,PLOT=TRUE,yoffset=-60000,xoffset = 50)
gam_lookup(DATASET = "HCT",DATE_HCT = 425,DATE_DEATH = NA,DATE_LDE = 1825,PLOT=TRUE,yoffset=-60000,xoffset = 50)
gam_lookup(DATASET = "HCT",DATE_HCT = 60,DATE_DEATH = 960,DATE_LDE = NA,PLOT=TRUE,yoffset=-60000,xoffset = 50)
gam_lookup(DATASET = "HCT",DATE_HCT = 60,DATE_DEATH = NA,DATE_LDE = 960,PLOT=TRUE,yoffset=-60000,xoffset = 50)
dev.off()


# CHEMO
chemo = "C:/ chemo.pdf"
pdf(file= chemo)
par(mfrow=c(3,2), cex=.5)
gam_lookup(DATASET = "CHEMO",DATE_HCT = NA,DATE_DEATH = 1825,DATE_LDE = NA,PLOT=TRUE,yoffset= -60000,xoffset=50)
gam_lookup(DATASET = "CHEMO",DATE_HCT = NA,DATE_DEATH = NA,DATE_LDE = 1825,PLOT=TRUE,yoffset= -60000,xoffset = 50)
gam_lookup(DATASET = "CHEMO",DATE_HCT = NA,DATE_DEATH = 960,DATE_LDE = NA,PLOT=TRUE,yoffset=-60000,xoffset = 50)
gam_lookup(DATASET = "CHEMO",DATE_HCT = NA,DATE_DEATH = NA,DATE_LDE = 960,PLOT=TRUE,yoffset=-60000,xoffset = 50)
dev.off()

#Get bootstrap: Standard deviation and 95% CI for HCT
res.mat <- matrix(data=NA,nrow = 1000,ncol = length(as.numeric(hct.store.eq[[1]]$coefficients)))
for (j in 1:1000){
  res.mat[j,] <- as.numeric(hct.store.eq[[j+1]]$coefficients)
}

ci.mat <- matrix(data=NA,nrow=length(as.numeric(hct.store.eq[[1]]$coefficients)),ncol=3)
for (k in 1:length(as.numeric(hct.store.eq[[1]]$coefficients))){
  ci.mat[k,] <- c(sd(res.mat[,k]),quantile(res.mat[,k],c(.025,.975)))
}

#Get bootstrap: Standard deviation and 95% CI for Chemo
chemores.mat <- matrix(data=NA,nrow = 1000,ncol = length(as.numeric(chemo.store.eq[[1]]$coefficients)))
for (j in 1:1000){
  chemores.mat[j,] <- as.numeric(chemo.store.eq[[j+1]]$coefficients)
}

chemoci.mat <- matrix(data=NA,nrow=length(as.numeric(chemo.store.eq[[1]]$coefficients)),ncol=3)
for (k in 1:length(as.numeric(chemo.store.eq[[1]]$coefficients))){
  chemoci.mat[k,] <- c(sd(chemores.mat[,k]),quantile(chemores.mat[,k],c(.025,.975)))
}

