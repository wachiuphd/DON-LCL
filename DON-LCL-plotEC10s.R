library(tidyverse)
library(ggplot2)
library(parallel)
library(rstan)

DON_dat <- data.frame(read.csv("DON-LCLdata_test.csv"))

### EC10 distribution ###
cells <- unique(DON_dat$Cell.line)
pec10 <- list()
load("DON_stanfit.RData")
set.seed(314159)
fitparms_df <- as.data.frame(rstan::extract(stan_fit))
ec10names <- paste("ec10",seq(1,length(cells)),sep=".")
ec10_df <- fitparms_df[,ec10names]
names(ec10_df) <- as.character(cells)
write.csv(t(apply(ec10_df,2,quantile,prob=c(0.025,0.5,0.975))),
          file="Table-DON_EC10_CellLine_summary.csv")
ec10_df.long <- melt(ec10_df)
pec10 <- ggplot(ec10_df.long) + 
  geom_boxplot(aes(x=reorder(variable,value,median),y=value),
               outlier.shape = NA) +
  coord_flip() + scale_y_log10() + ylab("EC10 (uM)") + xlab("Cell Line") +
  theme(axis.title = element_text(size = 15))
print(pec10)
ggsave("Figure-DON_EC10.pdf", plot=pec10,width=8,height=10)

### EC10 Population GM and GSD values
load("DON_stanfit.RData")
set.seed(314159)
fitparms_df <- as.data.frame(rstan::extract(stan_fit))
ec10.pop <- exp(fitparms_df[,c("m_x10","sd_x10")])
names(ec10.pop) <- c("EC10.GM","EC10.GSD")
write.csv(ec10.pop,file="DON-EC10_pop_samples.csv")
write.csv(t(apply(ec10.pop,2,quantile,prob=c(0.025,0.5,0.975))),
          "Table-EC10.pop_posteriors.csv")
ggplot(ec10.pop)+geom_histogram(aes(x=EC10.GM,y=..density..))+
  scale_x_log10()+annotation_logticks(side="b")+theme_bw()
ggsave("Figure-DON-EC10_pop.GM.pdf",height=3,width=5)
ggplot(ec10.pop)+geom_histogram(aes(x=EC10.GSD,y=..density..))+
  scale_x_log10()+annotation_logticks(side="b")+theme_bw()
ggsave("Figure-DON-EC10_pop.GSD.pdf",height=3,width=5)
prior.dens <- data.frame(xlog10=seq(0,1.5,0.01))
prior.dens$density.gsd <- dlnorm(prior.dens$xlog10,meanlog=log(0.221),sdlog=log(1.89))
prior.dens$x.gsd <- 10^prior.dens$xlog10

ggplot(ec10.pop)+geom_histogram(aes(x=EC10.GSD,y=..density..))+
  geom_line(aes(x.gsd,density.gsd),data=prior.dens)+
  geom_vline(xintercept=10^(0.5/qnorm(0.99)),color="red",linetype="dotted")+
  scale_x_log10(limits=c(1,10))+annotation_logticks(side="b")+theme_bw()
ggsave("Figure-DON-EC10_pop.GSD_compare_prior.pdf",height=3,width=5)

prior.dens$x.tdvf05 <- prior.dens$x.gsd^qnorm(0.95)
prior.dens$density.tdvf05 <- prior.dens$density.gsd/qnorm(0.95)
ec10.pop$TDVF05 <- ec10.pop$EC10.GSD^qnorm(0.95)
ggplot(ec10.pop)+geom_histogram(aes(x=TDVF05,y=..density..))+
  geom_line(aes(x.tdvf05,density.tdvf05),data=prior.dens)+
  geom_vline(xintercept=10^(0.5),color="red",linetype="dotted")+
  scale_x_log10(limits=c(1,30))+annotation_logticks(side="b")+theme_bw()
ggsave("Figure-DON-EC10_pop.TDVF05_compare_prior.pdf",height=3,width=5)
ec10.pop$TDVF01 <- ec10.pop$EC10.GSD^qnorm(0.99)
prior.dens$x.tdvf01 <- prior.dens$x.gsd^qnorm(0.99)
prior.dens$density.tdvf01 <- prior.dens$density.gsd/qnorm(0.99)
ggplot(ec10.pop)+geom_histogram(aes(x=TDVF01,y=..density..))+
  geom_line(aes(x.tdvf01,density.tdvf01),data=prior.dens)+
  geom_vline(xintercept=10^(0.5),color="red",linetype="dotted")+
  scale_x_log10(limits=c(1,30))+annotation_logticks(side="b")+theme_bw()
ggsave("Figure-DON-EC10_pop.TDVF01_compare_prior.pdf",height=3,width=5)  

write.csv(t(apply(ec10.pop,2,quantile,prob=c(0.025,0.5,0.975))),
          "Table-EC10.pop_posteriors.csv")

### EC10 deciles ### not used
# 
# decilefunc <- function(v) {
#   d <- as.integer(cut(as.numeric(v),c(quantile(v,probs=(0:9/10)),Inf),
#                       include.lowest=TRUE))
#   return(d)
# }
# 
# countfunc <- function(v) {
#   c <- hist(v,plot=FALSE,breaks=0.5+0:10)$counts
#   return(c)
# }
# 
# cells <- unique(DON_dat$Cell.line)
# pec10dec <- list()
# load("DON_stanfit.RData")
# fitparms_df <- as.data.frame(rstan::extract(stan_fit))[1:500,]
# ec10names <- paste("ec10",seq(1,length(cells)),sep=".")
# ec10_df <- fitparms_df[,ec10names]
# names(ec10_df) <- as.character(cells)
# 
# ec10_df_dec <- data.frame(t(apply(ec10_df,1,decilefunc)))
# names(ec10_df_dec) <- cells
# ec10_df_dec_counts <- data.frame(t(apply(ec10_df_dec,2,countfunc)))
# names(ec10_df_dec_counts) <- paste("decile",1:10,sep="")
# ec10_df_maxdec <- data.frame(CellLine=cells,
#                              decile_max=apply(ec10_df_dec_counts,
#                                               1,which.max),
#                              freq_max=apply(ec10_df_dec_counts,1,max)/500
# )
# ec10_df <- melt(ec10_df)
# ec10_summary <- aggregate(. ~ variable,data=ec10_df,median)
# names(ec10_summary) <- c("CellLine","EC10_50")
# ec10_summary$EC10_05 <- aggregate(. ~ variable,
#                                   data=ec10_df,quantile,prob=0.05)$value
# ec10_summary$EC10_95 <- aggregate(. ~ variable,
#                                   data=ec10_df,quantile,prob=0.95)$value
# ec10_summary <- merge(ec10_summary,ec10_df_maxdec)
# ec10_summary <- ec10_summary[order(ec10_summary$EC10_50),]
# write.csv(ec10_summary,file="DON_EC10_summary_deciles.csv",
#           row.names = FALSE)
