## This is the clean/final version of the R script doing TB simulation
rbinom <- Vectorize(rbinom,c("prob")) ## Vectorize the rbinom function for tb simulation

## hhsim function is the function for simulating tb transmission for each household
## Parameter explanations: 
## data: the simulated data for each household
## tpar: parameter values for the transmission model
## beta: the parameter that defines "transmission force" based on the household size, i.e., the term (1/n)^beta in the transmission model, see Cauchemez et al. (2004)
## apar: parameter values for the progression to disease model, right now there are only three: baseline, female and adult
## activep: the probability of progressing to active disease 
hhsim = function(data,tpar,beta,apar){
  n <- nrow(data)
  transmat <- mat.or.vec(n,n)
  datemat <- mat.or.vec(n,n)
  data$id <- 1:nrow(data)
  lambdac <- tpar[unique(data$mun)]*(data$adult*tpar['adult']+(1-data$adult))*(data$female*tpar['female']+(1-data$female))
  lambdaa <- apar['baseline']*(data$adult*apar['adult']+(1-data$adult))*(data$female*apar['female']+(1-data$female))
  p <- 1-exp(-lambdac)
  activep <- 1-exp(-lambdaa)
  notb <- rbinom(1,1,prob = p)
  if(n==1){
    actb <- ifelse(notb>0,rbinom(1,1,prob = activep),0)
    data$infect <- notb + actb
    if (data$infect == 2){
      data$date <- runif(1,0,10)
    }
  } else {
    if (any(notb> 0)) {
      actb <- numeric(n)
      actb[which(notb > 0)] <- as.vector(rbinom(1,1,prob = activep[which(notb>0)]))
      trans <- notb + actb
    } else {
      trans <- notb
    }
    date <- numeric(n)
    diag(transmat) <- trans
    if (any(trans == 2)){
      date[which(trans==2)] <- runif(length(which(trans==2)),0,10)
      diag(datemat) <- date
      t <- min(date[date>0]) ## The diagnosing time of the sender, in this case it is the initial case
      actind <- which(date==t)
      ind <- 1:nrow(data)
      sender <- tail(actind,n=1)
      receiver <- ind[-actind]
      while(length(actind) <= n){
        crowd <- unique(data$crowd)
        adult <- data[receiver,]$adult
        female <- data[receiver,]$female
        culture <- data[sender,]$culture
        smear <- data[sender,]$smear
        lambdah <- tpar['theta']*n^(beta)*(crowd*tpar['crowd1']+(1-crowd))*(adult*tpar['adult']+(1-adult))*(female*tpar['female']+(1-female))
        lambdaa <- apar['baseline']*(adult*apar['adult']+(1-adult))*(female*apar['female']+(1-female))
        p <- 1-exp(-lambdah)
        activep <- 1-exp(-lambdaa)
        notb <- rbinom(1,1,p)
        if (any(notb > 0)){
          actb <- numeric(length(p))
          actb[which(notb>0)] <- as.vector(rbinom(1,1,activep[which(notb>0)]))
          trans <- actb+notb
        } else {
          trans <- notb
        }
        date <- numeric(length(p))
        if (any(trans==2)){
          date[which(trans==2)] <- t+rgamma(sum(trans==2),1.46,0.657)
        } 
        transmat[sender,receiver] <- trans
        datemat[sender,receiver] <- date
        if (all(trans < 2) & all(datemat[,-actind]==0)){break}
        k <- datemat[,-actind]
        t <- min(k[k>0])
        sender <- which(datemat==t,arr.ind = TRUE)[,2]
        actind <- c(actind,sender)
        if (length(actind) == n) {break}
        receiver <- ind[-actind]
      }
    }
    data$infect <- apply(transmat,2,function(x) max(x,na.rm = TRUE))
    data$date <- apply(datemat,2,function(x) ifelse(any(x>0),min(x[x>0]),NA))
    data <- data%>%arrange(desc(infect),date)
    transmat <- transmat[data$id,data$id]
    transmat[lower.tri(transmat)] <- NA
    com_noactive <- diag(transmat)[which(data$infect!=2)]
    transmat[which(data$infect!=2),which(data$infect!=2)] <- NA
    diag(transmat)[which(data$infect!=2)] <- com_noactive
    data$edge <- apply(transmat,2,function(x) paste(as.character(x)[!is.na(as.character(x))],collapse = ""))
    data <- data[,-7]
  }
  return(data)
}

hhsim <- compiler::cmpfun(hhsim)

## tbsim is the function we will actually use
## it basically do three things: 
## First, it will generate simulated data for the entire population
## Second, it will call hhsim to simulate tb transmission for each household
## Third, it will do some postprocessing, which generate data for household contact study (using a sampler)
## The parameter prevalence specifies the probability of crowdedness, adult, female, smear and culture in the simulation
## The parameter munsize specifies the (stratified) sampling probabilities of all the communities 
tbsim = function(hhno,hhmin,hhmean,tpar,beta,apar,prevalence,munsize){
  library(tidyverse)
  mypar <- tpar
  mybeta <- beta
  myactivep <- apar
  names(mypar) <- c('theta',c('CARIACICA','SERRA','VITORIA','VILA VELHA'),'crowd1','adult','female')
  names(myactivep) <- c('baseline','adult','female')
  names(prevalence) <- c('crowd','adult','female')
  hhsize <- rpois(hhno, lambda = hhmean)
  l <- length(which(hhsize<hhmin))
  hhsize[hhsize<hhmin] <- sample(4:9,size = l,replace = TRUE)
  hhid <- rep(1:hhno,hhsize)
  crowd <- rep(sample(0:1,hhno,replace=TRUE,prob = c(1-prevalence['crowd'],prevalence['crowd'])),times=hhsize)
  mun <- rep(sample(c('CARIACICA','SERRA','VITORIA','VILA VELHA'),hhno,replace=TRUE,prob = munsize),times=hhsize)
  np <- sum(hhsize)
  female <- rbinom(np,1,prob = prevalence['female'])
  adult <- rbinom(np,1,prob = prevalence['adult'])
  date <- rep(NA,np)
  dat <- cbind.data.frame(hhid,crowd,mun,female,adult,date,stringsAsFactors=FALSE)
  names(dat) <- c('hhid','crowd','mun','female','adult','date')
  dat <- split(dat,dat$hhid)
  output <- lapply(dat,'hhsim',tpar=mypar,beta=mybeta,apar=myactivep)
  output <- do.call(rbind,output)
  ## result <- output%>%group_by(hhid)%>%count(infect)%>%spread(key=infect,value=n,fill=0)
  ## names(result) <- c('hhid','NoTB','LTB','ATB')
  ## result <- data.frame(result)
  hh1 <- output%>%group_by(hhid)%>%summarize(inf=sum(infect==2))%>%filter(inf>0)
  hh1 <- hh1$hhid
  output <- output%>%mutate(hhc=ifelse(hhid%in%hh1,1,0))
  return(output)
}

tbsim <- compiler::cmpfun(tbsim)


############################################################
## Posterior Prediction: Use simulation to obtain results ##
############################################################

#########################
## Parallel Processing ##
#########################

library(parallel)
## First, read the posterior estimates
w <- readRDS('est.rds')
P <- detectCores(logical = FALSE) ## P = 4
cl <- makeCluster(P)
clusterEvalQ(cl, {rbinom <- Vectorize(rbinom,c("prob"));library(tidyverse)})
clusterExport(cl,c("tbsim","hhsim","w"))
prediction <- function(tpar1){
  dat <- tbsim(hhno=10000,hhmin=4,hhmean = 6,tpar=tpar1,beta=0,apar=c(0.1,1,1),prevalence = c(0.4,0.7,0.5),munsize = rep(0.25,4))
  library(tidyverse)
  dat <- dat%>%group_by(hhid)%>%mutate(lag = as.numeric(date-min(date,na.rm=TRUE)))
  dat <- dat%>%ungroup(hhid)
  dat <- dat%>%arrange(hhid,desc(infect),lag)
  ## How to deal with two cases diagnosed at the same day?
  ## Take a look at the household whose hhid is 53
  ## I intentionally change the data for id:SE301009 from 2005-06-09 to 2005-06-10
  ## dat[dat$HHC=='SE301009','lag'] <- 1/365
  dat <- dat%>%group_by(hhid)%>%mutate(nopar1=ifelse(is.na(lag),sum(!is.na(lag))+1,ifelse(lag==0,1,rank(lag,na.last = NA))),hhsize=n())
  dat <- dat%>%mutate(nopar=ifelse(infect==0,0,nopar1)) ## non-TB cases have 0 parameter
  dat <- dat%>%mutate(nopar=ifelse(nopar==1,0,nopar)) ## Initial active TB cases have 0 parameter
  dat <- dat%>%ungroup(hhid)
  dat <- dat%>%mutate(y=ifelse(infect==0,0,1))
  dat <- dat%>%mutate(ca=ifelse(mun=='CARIACICA',1,0),se=ifelse(mun=='SERRA',1,0),vi=ifelse(mun=='VITORIA',1,0),vv=ifelse(mun=='VILA VELHA',1,0))
  dat <- dat%>%select(-c("mun","date"))
  ## How to process the edge into a path vector
  ## First of all, extracting the edge first
  s <- dat$edge
  t <- sum(dat$nopar1)
  out <- numeric(t)
  s <- paste(s,sep = "",collapse = "")
  for (i in 1:t){out[i] = as.numeric(substr(s,i,i))}
  out <- ifelse(out>0,1,0)
  ss <- cumsum(dat$nopar1)
  ## out is the vector we need
  ## Output 1: 5 statistics
  ## Probability of infection in household and community
  ## Probability of infection in household (pih)
  m <- which(dat$nopar>1)
  k=0
  p=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):(ss[i]-1)]
    k <- k+sum(out1)
    p <- p+length(out1)
  }
  pih <- k/p
  ## Probability of infection in community (pic, pic1, pic2, pic3, pic4)
  k=0
  for (i in 1:nrow(dat)){
    k <- k+out[ss[i]]
  }
  pic <- k/nrow(dat)
  rrhc <- pih/pic
  ## For each community
  k=0
  m <- which(dat$ca==1)
  for (i in m){
    k <- k+out[ss[i]]
  }
  pic1 <- k/length(m)
  m <- which(dat$nopar>1&dat$ca==1)
  k=0
  p=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):(ss[i]-1)]
    k <- k+sum(out1)
    p <- p+length(out1)
  }
  pih1 <- k/p
  rrhc1 <- pih1/pic1
  k=0
  m <- which(dat$se==1)
  for (i in m){
    k <- k+out[ss[i]]
  }
  pic2 <- k/length(m)
  m <- which(dat$nopar>1&dat$se==1)
  k=0
  p=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):(ss[i]-1)]
    k <- k+sum(out1)
    p <- p+length(out1)
  }
  pih2 <- k/p
  rrhc2 <- pih2/pic2
  k=0
  m <- which(dat$vi==1)
  for (i in m){
    k <- k+out[ss[i]]
  }
  pic3 <- k/length(m)
  m <- which(dat$nopar>1&dat$vi==1)
  k=0
  p=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):(ss[i]-1)]
    k <- k+sum(out1)
    p <- p+length(out1)
  }
  pih3 <- k/p
  rrhc3 <- pih3/pic3
  k=0
  m <- which(dat$vv==1)
  for (i in m){
    k <- k+out[ss[i]]
  }
  pic4 <- k/length(m)
  m <- which(dat$nopar>1&dat$vv==1)
  k=0
  p=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):(ss[i]-1)]
    k <- k+sum(out1)
    p <- p+length(out1)
  }
  pih4 <- k/p
  rrhc4 <- pih4/pic4
  ## Output 2: 3 statistics
  ## The relative susceptibility of adult (rsa)
  rsa <- (sum(dat$y[which(dat$adult==1)])/sum(dat$adult))/(sum(dat$y[which(dat$adult==0)])/sum(1-dat$adult))
  ## The relative susceptibility of female (rsf)
  rsf <- (sum(dat$y[which(dat$female==1)])/sum(dat$female))/(sum(dat$y[which(dat$female==0)])/sum(1-dat$female))
  ## The relative susceptibility of crowded family (rsc)
  rsc <- (sum(dat$y[which(dat$crowd==1)])/sum(dat$crowd))/(sum(dat$y[which(dat$crowd==0)])/sum(1-dat$crowd))
  ## Output 3: 10 statistics
  ## Still need to use the whole data
  ## The proportion of infection due to community transmission-Overall (prop_c)
  m <- which(dat$nopar>0)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)==1&out1[dat$nopar[i]]==1){k=k+1}
  }
  prop_c <- k/length(m) ## This is the proprotion we are looking for 
  ## The proportion of infection due to household transmission-Overall (prop_h)
  m <- which(dat$nopar>0)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)>0&out1[dat$nopar[i]]==0){k <- k+1}
  }
  prop_h <- k/length(m)
  ## For different communities
  ## The proportion of infection due to community transmission-1
  m <- which(dat$nopar>0&dat$ca==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)==1&out1[dat$nopar[i]]==1){k=k+1}
  }
  prop_c1 <- k/length(m) ## This is the proprotion we are looking for 
  ## The proportion of infection due to household transmission-1
  m <- which(dat$nopar>0&dat$ca==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)>0&out1[dat$nopar[i]]==0){k <- k+1}
  }
  prop_h1 <- k/length(m)
  ## The proportion of infection due to community transmission-2
  m <- which(dat$nopar>0&dat$se==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)==1&out1[dat$nopar[i]]==1){k=k+1}
  }
  prop_c2 <- k/length(m) ## This is the proprotion we are looking for 
  ## The proportion of infection due to household transmission-2
  m <- which(dat$nopar>0&dat$se==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)>0&out1[dat$nopar[i]]==0){k <- k+1}
  }
  prop_h2 <- k/length(m)
  ## The proportion of infection due to community transmission-3
  m <- which(dat$nopar>0&dat$vi==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)==1&out1[dat$nopar[i]]==1){k=k+1}
  }
  prop_c3 <- k/length(m) ## This is the proprotion we are looking for 
  ## The proportion of infection due to household transmission-3
  m <- which(dat$nopar>0&dat$vi==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)>0&out1[dat$nopar[i]]==0){k <- k+1}
  }
  prop_h3 <- k/length(m)
  ## The proportion of infection due to community transmission-4
  m <- which(dat$nopar>0&dat$vv==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)==1&out1[dat$nopar[i]]==1){k=k+1}
  }
  prop_c4 <- k/length(m) ## This is the proprotion we are looking for 
  ## The proportion of infection due to household transmission-4
  m <- which(dat$nopar>0&dat$vv==1)
  k=0
  for (i in m){
    out1 <- out[(ss[i-1]+1):ss[i]]
    if(sum(out1)>0&out1[dat$nopar[i]]==0){k <- k+1}
  }
  prop_h4 <- k/length(m)
  result <- c(rrhc,rrhc1,rrhc2,rrhc3,rrhc4,rsa,rsf,rsc,prop_c,prop_h,prop_c1,prop_h1,prop_c2,prop_h2,
              prop_c3,prop_h3,prop_c4,prop_h4)
  return(result)
}

out <- parApply(cl,w,1,prediction)
out <- t(out)
out <- data.frame(out)
colnames(out) <- c('rrhc','rrhc1','rrhc2','rrhc3','rrhc4','rsa','rsf','rsc','prop_c','prop_h','prop_c1','prop_h1','prop_c2','prop_h2',
                   'prop_c3','prop_h3','prop_c4','prop_h4')
saveRDS(out,'prediction.rds')
stopCluster(cl)










