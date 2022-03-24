library(tidyverse)
library(ggplot2)
library(parallel)
library(rstan)

drfunc <- function(x,parms) {
  yp <- parms$y0/ (1 + (x / parms$x0)^parms$n)
  return(yp)
}
parm_trans <- function(parms_sc) {
  y0 <- exp(parms_sc$m_y0 + parms_sc$sd_y0 * parms_sc$z_y0);
  x10 <- exp(parms_sc$m_x10 + parms_sc$sd_x10 * parms_sc$z_x10);
  n <- exp(parms_sc$m_n + parms_sc$sd_n * parms_sc$z_n);
  x0 <- x10 * 9.0^(1/n)
  return(data.frame(y0=y0,x10=x10,n=n,x0=x0))
}

DON_dat <- data.frame(read.csv("DON-LCLdata_test.csv"))
p <- ggplot(DON_dat)+
  geom_point(aes(x=conc,y=value))+
  scale_x_log10()+facet_wrap(~Cell.line,ncol=10)
print(p)

###plot fits###
xplot<-10^(seq(-10,35)/10)
cells <- unique(DON_dat$Cell.line)
pfit <- list()
source("DONstan_dat.R")
data_df <- data.frame(cell=cell,x=x,ys=ys,CellLine=cells[cell])

# # Check matches original data
# ggplot(onechemdat)+
#   geom_point(aes(x=conc,y=value/scale_factor,color=Chem))+
#   geom_point(data=data_df,aes(x=x,y=ys),size=0.1) +
#   scale_x_log10()+facet_wrap(~CellLine,ncol=10)

pred_df <- data.frame()
load("DON_stanfit.RData")
set.seed(314159)
fitparms_df <- as.data.frame(rstan::extract(stan_fit))[1:500,]
for (cellnum in 1:length(unique(cell))) {
  cat(cellnum,"...",sep="")
  parms_sc <- fitparms_df[,c("m_y0","m_x10","m_n","sd_y0","sd_x10","sd_n",
                             paste(c("z_y0","z_x10","z_n"),cellnum,sep="."))]
  parms <- parm_trans(parms_sc)
  cell_dat <- subset(data_df,CellLine==cells[cellnum])
  dfiter <- data.frame()
  for (j in 1:nrow(parms)) {
    dfiter <-rbind(dfiter,
                   data.frame(iter=j,x=xplot,yp=drfunc(xplot,parms[j,])))
  }
  pred_sum <- data.frame(CellLine=cells[cellnum],
                         x=xplot,
                         y50=aggregate(.~x,dfiter[,-1],median)$yp,
                         y025=aggregate(.~x,dfiter[,-1],quantile,0.025)$yp,
                         y975=aggregate(.~x,dfiter[,-1],quantile,0.975)$yp)
  pred_df <- rbind(pred_df,pred_sum)
}

cat("\n")
pfit <- ggplot(pred_df) + 
  scale_x_log10() +
  geom_ribbon(data=pred_df,aes(x=x,ymin=y025,ymax=y975),
              fill="grey",color="grey")+
  geom_line(data=pred_df,aes(x=x,y=y50),color="red") +
  geom_point(aes(x=x,y=ys),data=data_df,size=0.2) + 
  theme(axis.text.x = element_text(angle = 90), 
        axis.title = element_text(size = 15))+
  facet_wrap(~CellLine,ncol=11)+
  xlab("Concentration (uM)")+
  ylab("Viability")
print(pfit)
ggsave("Figure-DON_fit.pdf", plot=pfit,width=8, height=10)
