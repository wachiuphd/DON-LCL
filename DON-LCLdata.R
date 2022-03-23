###reformat the DON-LCLdata.csv###
library(tidyr)
library(dplyr)
library(readr)

data <- read.csv("DON-LCLdata.csv", fileEncoding = "UTF-8-BOM")
data_reformat <- pivot_longer(
  data, cols= dose_0.1_r1:dose_1000_r4, 
  names_to = "variable", values_to = "value")

library(stringi)
DON <- mutate(data_reformat, "conc" = case_when(
  stri_detect_fixed(variable, "dose_0.1_") ~ 0.1,
  stri_detect_fixed(variable, "dose_1_") ~ 1,
  stri_detect_fixed(variable, "dose_3_") ~ 3,
  stri_detect_fixed(variable, "dose_10_") ~ 10,
  stri_detect_fixed(variable, "dose_30_") ~ 30,
  stri_detect_fixed(variable, "dose_100_") ~ 100,
  stri_detect_fixed(variable, "dose_1000_") ~ 1000,
))

write.csv(DON,file="DON-LCLdata_test.csv",row.names = F)

###Prepare data for stan fitting###
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(parallel)
library(rstan)
library(MASS)
library(Hmisc)
library(lattice)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

DON_dat <- data.frame(read.csv("DON-LCLdata_test.csv"))
DON_dat$Cell.line <- factor(DON_dat$Cell.line)
ggplot(DON_dat)+  geom_point(aes(x=conc,y=value))+
  scale_x_log10()+facet_wrap(~Cell.line,ncol=10)

quants <- c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
xplot<-10^(seq(-10,35)/10)
scale_factor <- 100

DONstan_dat <- list(
    scale_factor = scale_factor,
    Ni = length(unique(DON_dat$Cell.line)),
    Nj = nrow(DON_dat),
    x = DON_dat$conc,
    ys = DON_dat$value/scale_factor,
    cell = as.numeric(DON_dat$Cell.line),
    quants = quants,
    Nquants = length(quants)
  )
 DON_fname <- "DONstan_dat.R"
 with(DONstan_dat, {
    stan_rdump(names(DONstan_dat),file=DON_fname)
 })
  
###rstan fitting###
 seed = 314159
 iter = 1000
 # fileprefix <- "DON"
 
 DONstan_dat <- list(
   scale_factor = scale_factor,
   Ni = length(unique(DON_dat$Cell.line)),
   Nj = nrow(DON_dat),
   x = DON_dat$conc,
   ys = DON_dat$value/scale_factor,
   cell = as.numeric(DON_dat$Cell.line),
   quants = quants,
   Nquants = length(quants)
 )  
 source("DONstan_dat.R")
 time.start<-proc.time()
 stan_fit <- stan(file="conc_resp_zero_ec10.stan",
                  data = DONstan_dat,
                  seed=seed,iter=iter, 
                  chains = 4,
                  sample_file="DONstan_dat_samples")
 time.end<-proc.time()
 print(time.end-time.start)
 save(stan_fit,file="DON_stanfit.RData")

 rhat.list <- stan_rhat(stan_fit)
 pdf(file="DON.pdf")
 plot(rhat.list)
 dev.off()

###Results###
   knitr::opts_chunk$set(echo = TRUE)
   library(tidyverse)
   library(readxl)
   library(reshape2)
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
ggsave("DON_LCL_data.pdf",plot=p,width=8,height=10)

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
    xlab("Concentration (mmol/L)")+
    ylab("Viability")
  print(pfit)
  ggsave("DON_fit.pdf", plot=pfit,width=8, height=10)
  ggsave("DON_fit.png", plot=pfit,width=8, height=10)

### EC10 distribution ###
  cells <- unique(DON_dat$Cell.line)
  pec10 <- list()
  load("DON_stanfit.RData")
  fitparms_df <- as.data.frame(rstan::extract(stan_fit))[1:500,]
  ec10names <- paste("ec10",seq(1,length(cells)),sep=".")
  ec10_df <- fitparms_df[,ec10names]
  names(ec10_df) <- as.character(cells)
  ec10_df <- melt(ec10_df)
  pec10 <- ggplot(ec10_df) + 
    geom_boxplot(aes(x=reorder(variable,value,median),y=value),
                   outlier.shape = NA) +
      coord_flip() + scale_y_log10() + ylab("EC10") + xlab("Cell Line") +
      theme(axis.title = element_text(size = 15))
  print(pec10)
  ggsave("DON_EC10.pdf", plot=pec10,width=8,height=10)
  ggsave("DON_EC10.png", plot=pec10,width=8,height=10)
  ec10_summary <- aggregate(. ~ variable,data=ec10_df,median)
  names(ec10_summary) <- c("CellLine","EC10_50")
  ec10_summary$EC10_05 <- aggregate(. ~ variable,
                                    data=ec10_df,quantile,prob=0.05)$value
  ec10_summary$EC10_95 <- aggregate(. ~ variable,
                                    data=ec10_df,quantile,prob=0.95)$value
  ec10_summary <- ec10_summary[order(ec10_summary$EC10_50),]
  write.csv(ec10_summary,file="DON_EC10_summary.csv",
            row.names = FALSE)

### EC10 Population GM and GSD values
  ec10.gsd <- exp(as.data.frame(extract(stan_fit,"sd_x10")))
  ec10.gm <- exp(as.data.frame(extract(stan_fit,"m_x10")))
  ec10.pop <- cbind(ec10.gm,ec10.gsd)
  names(ec10.pop) <- c("EC10.GM","EC10.GSD")
  write.csv(ec10.pop,file="EC10_pop.csv")
  write.csv(t(apply(ec10.pop,2,quantile,prob=c(0.025,0.5,0.975))),
            "Table-EC10.pop_posteriors.csv")
  ggplot(ec10.pop)+geom_histogram(aes(x=EC10.GM,y=..density..))+
    scale_x_log10()+annotation_logticks(side="b")+theme_bw()
  ggsave("EC10_pop.GM.pdf",height=3,width=5)
  ggplot(ec10.pop)+geom_histogram(aes(x=EC10.GSD,y=..density..))+
    scale_x_log10()+annotation_logticks(side="b")+theme_bw()
  ggsave("EC10_pop.GSD.pdf",height=3,width=5)
  
### EC10 deciles ###
    decilefunc <- function(v) {
      d <- as.integer(cut(as.numeric(v),c(quantile(v,probs=(0:9/10)),Inf),
                          include.lowest=TRUE))
      return(d)
    }
    
    countfunc <- function(v) {
      c <- hist(v,plot=FALSE,breaks=0.5+0:10)$counts
      return(c)
    }
    
    cells <- unique(DON_dat$Cell.line)
    pec10dec <- list()
    load("DON_stanfit.RData")
    fitparms_df <- as.data.frame(rstan::extract(stan_fit))[1:500,]
    ec10names <- paste("ec10",seq(1,length(cells)),sep=".")
    ec10_df <- fitparms_df[,ec10names]
    names(ec10_df) <- as.character(cells)
      
    ec10_df_dec <- data.frame(t(apply(ec10_df,1,decilefunc)))
    names(ec10_df_dec) <- cells
    ec10_df_dec_counts <- data.frame(t(apply(ec10_df_dec,2,countfunc)))
    names(ec10_df_dec_counts) <- paste("decile",1:10,sep="")
    ec10_df_maxdec <- data.frame(CellLine=cells,
                                 decile_max=apply(ec10_df_dec_counts,
                                                  1,which.max),
                                 freq_max=apply(ec10_df_dec_counts,1,max)/500
    )
    ec10_df <- melt(ec10_df)
    ec10_summary <- aggregate(. ~ variable,data=ec10_df,median)
    names(ec10_summary) <- c("CellLine","EC10_50")
    ec10_summary$EC10_05 <- aggregate(. ~ variable,
                                      data=ec10_df,quantile,prob=0.05)$value
    ec10_summary$EC10_95 <- aggregate(. ~ variable,
                                      data=ec10_df,quantile,prob=0.95)$value
    ec10_summary <- merge(ec10_summary,ec10_df_maxdec)
    ec10_summary <- ec10_summary[order(ec10_summary$EC10_50),]
    write.csv(ec10_summary,file="DON_EC10_summary_deciles.csv",
              row.names = FALSE)
###extract TDVF values###
    load("DON_stanfit.Rdata")
    ec10.med.new <- extract(stan_fit,"ec10_quants[6]")
    tdvf.01.new <- extract(stan_fit,"ec10_quant_ratios[1]")
    tdvf.05.new <- extract(stan_fit,"ec10_quant_ratios[3]")
    ec10.gm.new <- extract(stan_fit,"m_x10")
    ec10.gsd.new <- extract(stan_fit,"sd_x10")

    ec10.med.new <- as.data.frame(ec10.med.new)
    tdvf.01.new <- as.data.frame(tdvf.01.new)
    tdvf.05.new <- as.data.frame(tdvf.05.new)
    ec10.gm.new <- exp(as.data.frame(ec10.gm.new))
    ec10.gsd.new <- exp(as.data.frame(ec10.gsd.new))
    tdvf.01.alt <- ec10.gsd.new^qnorm(0.99)
    tdvf.05.alt <- ec10.gsd.new^qnorm(0.95)
    
    colnames(tdvf.01.new) <- "TDVF.01"
    colnames(tdvf.05.new) <- "TDVF.05"
    colnames(tdvf.01.alt) <- "TDVF.01"
    colnames(tdvf.05.alt) <- "TDVF.05"
  #  ec10.med.new <- melt(ec10.med.new)
  #  ec10.med.new$group <- "EC10.Med"
  #  tdvf.01.new <- melt(tdvf.01.new)
  #  tdvf.01.new$group <- "TDVF.01"
  #  tdvf.05.new <- melt(tdvf.05.new)
  #  tdvf.05.new$group <- "TDVF.05"
    
    # ggplot(tdvf.01.new) + geom_boxplot(aes(x=TDVF.01),outlier.shape = NA) +
    #   scale_x_log10() + xlab("TDVF01") + 
    #   annotation_logticks(sides="b") + coord_cartesian() + theme_classic()
    # ggsave("Replicates_boxplot_TDVF01.pdf")
    # ggsave("Replicates_boxplot_TDVF01.png")
    # 
    # ggplot(tdvf.05.new) + geom_boxplot(aes(x=TDVF.05),outlier.shape = NA) +
    #   scale_x_log10() + xlab("TDVF05") + 
    #   annotation_logticks(sides="b") + coord_cartesian() + theme_classic()
    # ggsave("Replicates_boxplot_TDVF05.pdf")
    # ggsave("Replicates_boxplot_TDVF05.png")
    
    ggplot(tdvf.01.new) + geom_histogram(aes(x=TDVF.01)) +
     scale_x_log10() + xlab("TDVF01") + ylab("DON") +
     annotation_logticks(sides="b") + coord_cartesian() + theme_classic()
    ggsave("Replicates_histogram_TDVF01.pdf")
    
    ggplot(tdvf.01.alt) + geom_histogram(aes(x=TDVF.01)) +
      scale_x_log10() + xlab("TDVF01") + ylab("DON") +
      annotation_logticks(sides="b") + coord_cartesian() + theme_classic()
    ggsave("Replicates_histogram_TDVF01.alt.pdf")
    
    ggplot(tdvf.05.new) + 
      geom_histogram(aes(x=TDVF.05), color="black", fill="white") +
      scale_x_log10(limits = c(1, 50), 
                    breaks=c(1, 3.16, 10, 20, 30, 40, 50), 
                    labels=c(1, 3.16, 10, 20, 30, 40, 50)) +
      scale_y_continuous(expand = c(0, 0))+
      xlab("TDVF05") + 
      annotation_logticks(sides="b") + coord_cartesian() + 
      theme_classic()+
      theme(axis.text.x = element_text(size = 15), 
            axis.text.y = element_text(size = 15), 
            axis.title = element_text(size = 20), 
            axis.title.y=element_blank())+
      geom_vline(xintercept=3.16, color="red", linetype="dotted")
    ggsave("Replicates_histogram_TDVF05.pdf", width=10, height=8)

    ggplot(tdvf.05.alt) + 
      geom_histogram(aes(x=TDVF.05), color="black", fill="white") +
      scale_x_log10(limits = c(1, 50), 
                    breaks=c(1, 3.16, 10, 20, 30, 40, 50), 
                    labels=c(1, 3.16, 10, 20, 30, 40, 50)) +
      scale_y_continuous(expand = c(0, 0))+
      xlab("TDVF05") + 
      annotation_logticks(sides="b") + coord_cartesian() + 
      theme_classic()+
      theme(axis.text.x = element_text(size = 15), 
            axis.text.y = element_text(size = 15), 
            axis.title = element_text(size = 20), 
            axis.title.y=element_blank())+
      geom_vline(xintercept=3.16, color="red", linetype="dotted")
    ggsave("Replicates_histogram_TDVF05.alt.pdf", width=10, height=8)
