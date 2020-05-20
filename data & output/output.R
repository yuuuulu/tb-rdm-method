## Read the DIC output and posterior sample
## The first 5 columns are weights favoring community transmission and the last 5 columns are weights favoring household transmission
out1 <- readRDS('out1.rds')
out2 <- readRDS('out2.rds')
out3 <- readRDS('out3.rds')
out4 <- readRDS('out4.rds')
## Need another function to quickly identify the optimal weights
## The input should be a vector (of DIC)
## Type = 1 is for extra-household transmission > household transmission
## Type = 0 is for household transmission > extra-household transmission
## This version is for knowing the direction
opweight <- function(dic,type){
  if (type == 1){
    s <- min(dic[1:5])
    k <- which(dic==s)
    ind <- which(dic<(s+4.5))
    final = ind[ind<6&ind>=k]
  } else if (type == 0) {
    s <- min(dic[6:10])
    k <- which(dic==s)
    ind <- which(dic<(s+4.5))
    final = ind[ind>5&ind<=k]
  }
  return(final)
}

out1 <- out1[,-6]
out2 <- out2[,-6]
out3 <- out3[,-6]
out4 <- out4[,-6]
bw1 <- apply(out1,1,opweight,type = 1) ## This should be the choice of weights based on DIC
bw2 <- apply(out2,1,opweight,type = 1)
bw3 <- apply(out3,1,opweight,type = 0)
bw4 <- apply(out4,1,opweight,type = 0)


## Get the estimates based on optimal weighting scheme
west <- function(sce_no,data,w){
  n <- length(w)
  file <- paste0("s",sce_no,"_est_",data,".rds")
  est <- readRDS(file)
  for (i in 1:n){
    if (i == 1){
      output <- est[[w[i]]]
    } else {
      output <- rbind(output,est[[w[i]]])
    }
  }
  out <- apply(output,2,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
  return(out)
}

west(1,1,bw1[[1]])

## wrap it up
wci <- function(sce_no,w){
  output <- vector("list",25)
  for (i in 1:25){
    output[[i]] <- west(sce_no,i,w[[i]])
  }
  return(output)
}
wci<- compiler::cmpfun(wci)
wci(1,bw1)

## DIC plot
par(mfrow=c(2,2))

dic <- apply(out1,2,mean)
plot(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic,type='l',ylab = 'DIC', xlab = 'The average weight of household transmission relative to extra-household transmission',main = 'Scenario 1: Predominant Extra-Household')
abline(h=min(dic)+5,lty=2,col='red')
points(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic)

dic <- apply(out2,2,mean)
plot(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic,type='l',ylab = 'DIC', xlab = 'The average weight of household transmission relative to extra-household transmission',main = 'Scenario 2: Stronger Extra-Household')
abline(h=min(dic)+5,lty=2,col='red')
points(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic)

dic <- apply(out3,2,mean)
plot(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic,type='l',ylab = 'DIC', xlab = 'The average weight of household transmission relative to extra-household transmission',main = 'Scenario 3: Stronger Household')
abline(h=min(dic)+5,lty=2,col='red')
points(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic)

dic <- apply(out4,2,mean)
plot(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic,type='l',ylab = 'DIC', xlab = 'The average weight of household transmission relative to extra-household transmission',main = 'Scenario 4: Predominant Household')
abline(h=min(dic)+5,lty=2,col='red')
points(c(0.78,0.82,0.86,0.91,0.95,1,1.02,1.04,1.06,1.08,1.1),dic)

par(mfrow=c(1,1))

## Compare the HHC and whole data first
library(tidyverse)
## Using function: s is the scenario number, best weight is the weight order, weightdic is the choice of weights based on DIC
plotdata <- function(s){
  filename <- paste0('s',s,'_whole.rds')
  sw <- readRDS(filename)
  lb <- matrix(unlist(lapply(sw,function(x) apply(x,2,function(y) quantile(y,probs = 0.025)))),nrow = 25,byrow = T)
  lb <- apply(lb,2,mean)
  ub <- matrix(unlist(lapply(sw,function(x) apply(x,2,function(y) quantile(y,probs = 0.975)))),nrow = 25,byrow = T)
  ub <- apply(ub,2,mean)
  md <- matrix(unlist(lapply(sw,function(x) apply(x,2,function(y) quantile(y,probs = 0.5)))),nrow = 25,byrow = T)
  md <- apply(md,2,mean)
  sw_plot <- rbind(lb,md,ub)
  sw_plot <- data.frame(t(sw_plot))
  sw_plot$type <- 'Whole'
  ## The HHC data
  filename <- paste0('s',s,'_hhc.rds')
  sh <- readRDS(filename)
  lb <- matrix(unlist(lapply(sh,function(x) apply(x,2,function(y) quantile(y,probs = 0.025)))),nrow = 25,byrow = T)
  lb <- apply(lb,2,mean)
  ub <- matrix(unlist(lapply(sh,function(x) apply(x,2,function(y) quantile(y,probs = 0.975)))),nrow = 25,byrow = T)
  ub <- apply(ub,2,mean)
  md <- matrix(unlist(lapply(sh,function(x) apply(x,2,function(y) quantile(y,probs = 0.5)))),nrow = 25,byrow = T)
  md <- apply(md,2,mean)
  sh_plot <- rbind(lb,md,ub)
  sh_plot <- data.frame(t(sh_plot))
  sh_plot$type <- 'HHC'
  ## The HHC data with weights 
  filename <- paste0('s',s,'_hw.rds')
  shw <- readRDS(filename)
  lb <- matrix(unlist(lapply(shw,function(x) x[1,])),nrow = 25,byrow = T)
  lb <- apply(lb,2,mean)
  ub <- matrix(unlist(lapply(shw,function(x) x[3,])),nrow = 25,byrow = T)
  ub <- apply(ub,2,mean)
  md <- matrix(unlist(lapply(shw,function(x) x[2,])),nrow = 25,byrow = T)
  md <- apply(md,2,mean)
  shw_plot <- rbind(lb,md,ub)
  shw_plot <- data.frame(t(shw_plot))
  shw_plot$type <- 'HHC with weights'
  s_plot <- rbind(sw_plot,sh_plot,shw_plot)
  return(s_plot)
}

s1_plot <- plotdata(1)
s1_plot$scenario <- 'Predominant Extra-Household Transmission'
tv1 <- c(0.05,0.99,1.07,0.83,1.26,1.2,1.2,1)

s2_plot <- plotdata(2)
s2_plot$scenario <- 'Stronger Extra-Household Transmission'
tv2 <- c(0.21,0.83,0.91,0.67,1.08,1.2,1.2,1)

s3_plot <- plotdata(3)
s3_plot$scenario <- 'Stronger Household Transmission'
tv3 <-  c(0.8,0.3,0.38,0.2,0.5,1.2,1.2,1) 

s4_plot <- plotdata(4)
s4_plot$scenario <- 'Predominant Household Transmission'
tv4 <- c(1.2,0.1,0.15,0.05,0.2,1.2,1.2,1)

s1p1 <- s1_plot[c(1:5,9:13,17:21),]
s1p1$parameter <- rep(c('Household','Community1','Community2','Community3','Community4'),3)
s1p2 <- s1_plot[c(6:8,14:16,22:24),]
s1p2$parameter <- rep(c('Crowd','Adult','Female'),3)
s2p1 <- s2_plot[c(1:5,9:13,17:21),]
s2p1$parameter <- rep(c('Household','Community1','Community2','Community3','Community4'),3)
s2p2 <- s2_plot[c(6:8,14:16,22:24),]
s2p2$parameter <- rep(c('Crowd','Adult','Female'),3)
s3p1 <- s3_plot[c(1:5,9:13,17:21),]
s3p1$parameter <- rep(c('Household','Community1','Community2','Community3','Community4'),3)
s3p2 <- s3_plot[c(6:8,14:16,22:24),]
s3p2$parameter <- rep(c('Crowd','Adult','Female'),3)
s4p1 <- s4_plot[c(1:5,9:13,17:21),]
s4p1$parameter <- rep(c('Household','Community1','Community2','Community3','Community4'),3)
s4p2 <- s4_plot[c(6:8,14:16,22:24),]
s4p2$parameter <- rep(c('Crowd','Adult','Female'),3)


dat1 <- rbind(s1p1,s2p1,s3p1,s4p1)
dat2 <- rbind(s1p2,s2p2,s3p2,s4p2)
dat1$truevalue <-  c(rep(tv1[1:5],3),rep(tv2[1:5],3),rep(tv3[1:5],3),rep(tv4[1:5],3))
dat2$truevalue <-  c(rep(tv1[6:8],3),rep(tv2[6:8],3),rep(tv3[6:8],3),rep(tv4[6:8],3))


## Figure 1
ggplot(data = dat1)+geom_linerange(mapping = aes(x=parameter,ymin=lb,ymax=ub,group=type),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=parameter,y=md,shape = type),position = position_dodge(width = 0.4),size = 1.6)+geom_point(mapping = aes(x=parameter,y=truevalue), position = position_dodge(width = 0.4),size = 1.6, colour = 'red', shape = 4)+
  labs(x = "Parameter", y = "Estimates",shape = "Mean of Median" )+facet_wrap(~scenario, nrow = 2)


## Figure 2
ggplot(data = dat2)+geom_linerange(mapping = aes(x=parameter,ymin=lb,ymax=ub,group=type),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=parameter,y=md,shape = type),position = position_dodge(width = 0.4),size = 1.6)+geom_point(mapping = aes(x=parameter,y=truevalue), position = position_dodge(width = 0.4),size = 1.6, colour = 'red', shape = 4)+
  labs(x = "Parameter", y = "Estimates",shape = "Mean of Median" )+facet_wrap(~scenario, nrow = 2)


## Making Tables

coverage <- function(data,tv){
  lb <- data[1,] 
  ub <- data[3,] 
  x <- numeric(length(tv))
  for (i in 1:length(tv)){
    x[i] = ifelse(tv[i]<=ub[i] & tv[i] >= lb[i], 1, 0)
  }
  return(x)
}

table <- function(s,tv){
  filename <- paste0('s',s,'_whole.rds')
  sw <- readRDS(filename)
  filename <- paste0('s',s,'_hhc.rds')
  sh <- readRDS(filename)
  output <- mat.or.vec(3,8)
  sw <- lapply(sw,function(x) apply(x,2,function(y) quantile(y,probs = c(0.025,0.5,0.975))))
  o <- matrix(unlist(lapply(sw,'coverage',tv = tv)),nrow = 25, byrow = T)
  output[1,] <- apply(o,2,sum)
  sh <- lapply(sh,function(x) apply(x,2,function(y) quantile(y,probs = c(0.025,0.5,0.975))))
  o <- matrix(unlist(lapply(sh,'coverage',tv = tv)),nrow = 25, byrow = T)
  output[2,] <- apply(o,2,sum)
  filename <- paste0('s',s,'_hw.rds')
  shw <- readRDS(filename)
  o <- matrix(unlist(lapply(shw,'coverage',tv = tv)),nrow = 25, byrow = T)
  output[3,] <- apply(o,2,sum)
  rownames(output) <- c('Whole','HHC','HHC_W')
  return(output)
}

table(1,tv1)
table(2,tv2)
table(3,tv3)
table(4,tv4)

################################################################### 
## Simulation with powerful predictors and reparameterized model ##
###################################################################

out <- readRDS('powersim.rds')
ci <- lapply(out,function(x) apply(x,2,function(y) quantile(y,probs = c(0.025,0.5,0.975))))
coverage <- function(data,tv){
  lb <- data[1,] 
  ub <- data[3,] 
  x <- numeric(length(tv))
  for (i in 1:length(tv)){
    x[i] = ifelse(tv[i]<=ub[i] & tv[i] >= lb[i], 1, 0)
  }
  return(x)
}
## Make a table
o <- matrix(unlist(lapply(ci,'coverage',tv =c(0.03,1.2,1,1.5,5.5,2,2.8,3.5,3,2.5))),nrow = 1000, byrow = T)
apply(o,2,sum)
md <- matrix(unlist(lapply(ci,function(x) x[2,])),nrow=1000,byrow=T)
apply(md,2,mean)

## Make a figure
plotdata <- function(dat,size){
  lb <- matrix(unlist(lapply(dat,function(x) apply(x,2,function(y) quantile(y,probs = 0.025)))),nrow = size,byrow = T)
  lb <- apply(lb,2,mean)
  ub <- matrix(unlist(lapply(dat,function(x) apply(x,2,function(y) quantile(y,probs = 0.975)))),nrow = size,byrow = T)
  ub <- apply(ub,2,mean)
  md <- matrix(unlist(lapply(dat,function(x) apply(x,2,function(y) quantile(y,probs = 0.5)))),nrow = size,byrow = T)
  md <- apply(md,2,mean)
  plotd <- rbind(lb,md,ub)
  plotd <- data.frame(t(plotd))
  return(plotd)
}

d <- plotdata(out,1000)
tv <- c(0.03,1.2,1,1.5,5.5,2,2.8,3.5,3,2.5)
d$truevalue <- tv
d$parameter <- c('Baseline','Adult','Female','Crowd','Sleep','Time','Severity','Burden','SES','Cough')

ggplot(data = d)+geom_linerange(mapping = aes(x=parameter,ymin=lb,ymax=ub))+
  geom_point(mapping = aes(x=parameter,y=md),size = 1.8)+geom_point(mapping = aes(x=parameter,y=truevalue), position = position_dodge(width = 0.3),size = 1.8, colour = 'red', shape = 4)+
  labs(x = "Parameter", y = "Estimates" )+geom_hline(yintercept = 1,linetype = "dashed")



#########################
## Brazil Data Example ##
#########################
library(tidyverse)
library(wselect)
data <- read_csv("mydat.csv")
data[904,]$adult=1
data <- data%>%group_by(hhid)%>%mutate(lag = as.numeric((date-min(date,na.rm=TRUE))/365))
data <- data%>%ungroup(hhid)
data <- data%>%arrange(hhid,desc(infect),lag)
data <- data%>%group_by(hhid)%>%mutate(nopar1=ifelse(is.na(lag),sum(!is.na(lag))+1,ifelse(lag==0,1,rank(lag,na.last = NA))),hhsize=n())
data <- data%>%group_by(hhid)%>%mutate(nopar=ifelse(is.na(lag),sum(!is.na(lag))+1,ifelse(lag==0,1,rank(lag,na.last = NA))))
data <- data%>%mutate(nopar=ifelse(infect==0,0,nopar)) ## non-TB cases have 0 parameter
data <- data%>%mutate(nopar=ifelse(infect==2&nopar==1,0,nopar)) ## Initial active TB cases have 0 parameter
data <- data%>%ungroup(hhid)
data <- data%>%mutate(ca=ifelse(mun=='CARIACICA',1,0),se=ifelse(mun=='SERRA',1,0),vi=ifelse(mun=='VITORIA',1,0),vv=ifelse(mun=='VILA VELHA',1,0))
data <- data%>%select(-c("HHC","mun","date","hhsize"))
hhsize <- as.numeric(table(data$hhid))
hhsize <- c(1,hhsize[-160])
index <- cumsum(hhsize)
data1 <- data[-index,]
ts <- data1$infect
ts <- ifelse(ts>0,1,0)
out <- readRDS('weight_DIC.rds')
w <- readRDS('weight_est.rds')
u <- readRDS('unweight_est.rds')
mw <- lapply(w,function(x) apply(x,2,mean))
mw <- do.call(rbind,mw)
mu <- apply(u,2,mean)
## dic_mw <- apply(mw,1,adcheck,data=data,iter=2000,tv = ts)
## dic_mu <- adcheck(data=data,iter=2000,tv = ts,tpar=mu)
dic_mw <- c(-569.2417,-569.7481,-569.5295,-572.5498,-575.1341,-580.5582,-583.1833,-587.8702,-586.6880,-589.2513)
dic_mu = -577.9769
result <- c(dic_mw[1:5],dic_mu,dic_mw[6:10])
result <- -2*result ## Now result is the DIC(parameter_mean_value)
out <- -2*out 
result1 <- apply(out,2,mean) ## result1 is mean(DIC(parameter_value))
dic <- 2*result1-result
plot(c(0.78,0.82,0.86,0.91,0.95,1,1.05,1.1,1.16,1.22,1.28),dic,type='l',ylab = 'DIC', xlab = 'The average weight of household transmission relative to extra-household transmission')
abline(h=min(dic)+5,lty=2,col='red')
points(0.82,dic[2],pch=0)
points(c(0.86,0.91),dic[c(3,4)],pch=2)
points(c(0.78,0.95,1,1.05,1.1,1.16,1.22,1.28),dic[-c(2,3,4)])
## Conclusion: we should choose -0.2,-0.15 and -0.1 as the fitting weighting schemes
est <- rbind(w[[2]],w[[3]],w[[4]])
apply(est,2,function(x) quantile(x, probs = c(0.025,0.5,0.975)))
## Create Brazil figure 1
lb <- apply(u,2,function(x) quantile(x, probs = 0.025))
md <- apply(u,2,function(x) quantile(x, probs = 0.5))
ub <- apply(u,2,function(x) quantile(x, probs = 0.975))
parameter <- c('Household','Community1','Community2','Community3','Community4','Crowd','Adult','Female')
estimator <- 'HHC'
dat1 <- cbind.data.frame(lb,md,ub,parameter,estimator)
lb <- apply(est,2,function(x) quantile(x, probs = 0.025))
md <- apply(est,2,function(x) quantile(x, probs = 0.5))
ub <- apply(est,2,function(x) quantile(x, probs = 0.975))
parameter <- c('Household','Community1','Community2','Community3','Community4','Crowd','Adult','Female')
estimator <- 'HHC with weights'
dat2 <- cbind.data.frame(lb,md,ub,parameter,estimator)
dat <- rbind(dat1,dat2)

dat1 <- dat[c(1:5,9:13),]
dat2 <- dat[c(6:8,14:16),]

ggplot(data = dat1)+geom_linerange(mapping = aes(x=parameter,ymin=lb,ymax=ub,group=estimator),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=parameter,y=md,shape = estimator),position = position_dodge(width = 0.4),size = 1.8)+
  labs(title = "Estimates based on Brazilian household contact study",
       x = "Parameter", y = "Estimates",shape = "Median" )

ggplot(data = dat2)+geom_linerange(mapping = aes(x=parameter,ymin=lb,ymax=ub,group=estimator),position = position_dodge(width = 0.4))+
  geom_point(mapping = aes(x=parameter,y=md,shape = estimator),position = position_dodge(width = 0.4),size = 1.8)+
  labs(title = "Estimates based on Brazilian household contact study",
       x = "Parameter", y = "Estimates",shape = "Median" )+geom_hline(yintercept = 1,linetype = "dashed")

## Model prediction
## Interpretable results for the Brazil example
library(tidyverse)
out <- readRDS('prediction.rds')
apply(out,2,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
## Violin plot
wsq <- out[,1:5]
wsq <- wsq%>%gather(1:5,key='type',value = 'number')
ggplot(data=wsq,aes(x=type,y=number))+geom_violin()+scale_y_continuous(breaks = seq(0,1.5,by=0.1) ,name = 'Relative risk of household transmission versus extra-household transmission')+ 
  scale_x_discrete(breaks = unique(wsq$type),labels=c('Overall','Community1','Community2','Community3','Community4'),name=NULL)
wsq <- out[,6:8]
wsq <- wsq%>%gather(1:3,key='type',value = 'number')
ggplot(data=wsq,aes(x=type,y=number))+geom_violin()+scale_y_continuous(breaks = seq(0.8,1.4,by=0.1),name = 'Relative Risk of TB infection')+ 
  scale_x_discrete(breaks = unique(wsq$type),labels=c('Adults vs. Children','Female vs. Male','Crowded vs. Non-crowded'),name=NULL)+geom_hline(yintercept = 1,linetype = "dashed")
wsq <- out[,11:18]
wsq <- wsq%>%gather(1:8,key='type',value = 'number')
ggplot(data=wsq,aes(x=type,y=number))+geom_violin()+scale_y_continuous(breaks = seq(0,1,by=0.1),name = 'The proportion of infection attributed to')+ 
  scale_x_discrete(breaks = unique(wsq$type),labels=c('Extra-Household1','Household1','Extra-Household2','Household2','Extra-Household3','Household3','Extra-Household4','Household4'),name=NULL)

