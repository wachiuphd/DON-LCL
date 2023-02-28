library(tidyr)
library(dplyr)
library(readr)
library(stringi)
library(ggplot2)
library(parallel)
library(rstan)
library(MASS)
library(Hmisc)
library(lattice)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

###reformat the DON-LCLdata.csv###
data <- read.csv("DON-LCLdata.csv", fileEncoding = "UTF-8-BOM")
data_reformat <- pivot_longer(
  data, cols= dose_0.1_r1:dose_1000_r4, 
  names_to = "variable", values_to = "value")

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

DON_dat <- data.frame(read.csv("DON-LCLdata_test.csv"))
DON_dat$Cell.line <- factor(DON_dat$Cell.line)
ggplot(DON_dat)+  geom_point(aes(x=conc,y=value))+
  scale_x_log10()+facet_wrap(~Cell.line,ncol=10)

quants <- c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
xplot<-10^(seq(-10,35)/10)
scale_factor <- 100 #Emax

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
 iter = 2000
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
 pdf(file="DON-Rhat.pdf")
 plot(rhat.list)
 dev.off()

