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
  output <- output%>%group_by(hhid)%>%mutate(lag = as.numeric(date-min(date,na.rm=TRUE)))
  output <- output%>%ungroup(hhid)
  output <- output%>%arrange(hhid,desc(infect),lag)
  output <- output%>%group_by(hhid)%>%mutate(nopar1=ifelse(is.na(lag),sum(!is.na(lag))+1,ifelse(lag==0,1,rank(lag,na.last = NA))),hhsize=n())
  output <- output%>%mutate(nopar=ifelse(hhc==0,0,nopar1))
  output <- output%>%mutate(nopar=ifelse(hhc==1&infect==0,0,nopar)) ## non-TB cases have 0 parameter
  output <- output%>%mutate(nopar=ifelse(hhc==1&infect==2&nopar==1,0,nopar)) ## Initial active TB cases have 0 parameter
  output <- output%>%ungroup(hhid)
  output <- output%>%mutate(ca=ifelse(mun=='CARIACICA',1,0),se=ifelse(mun=='SERRA',1,0),vi=ifelse(mun=='VITORIA',1,0),vv=ifelse(mun=='VILA VELHA',1,0))
  output <- output%>%select(-c("mun","date"))
  return(output)
}

tbsim <- compiler::cmpfun(tbsim)

library(parallel)
P <- detectCores(logical = FALSE) ## P = 4
cl <- makeCluster(P)
clusterEvalQ(cl, {rbinom <- Vectorize(rbinom,c("prob"));library(tidyverse)})
clusterExport(cl,c("tbsim","hhsim"))
datagen <- function(tpar1,hhno1){
  result <- tbsim(hhno=hhno1,hhmin=4,hhmean = 6.5,tpar=tpar1,beta=0,apar=c(0.1,1,1),prevalence = c(0.4,0.7,0.5),munsize = rep(0.25,4))
  return(result)
}
## Simulate ten datasets per scenario 
## Scenario 1: Extremely Predominant Community (Over 80%)
## Use the prediction function in simulation_prediction R code to select the true parameter values
## Using tpar = c(0.05,0.99,1.07,0.83,1.26,1.2,1.2,1)
tpar <- matrix(rep(c(0.05,0.99,1.07,0.83,1.26,1.2,1.2,1),25),nrow = 25, byrow = T)
out <- parApply(cl,tpar,1,datagen,hhno1=500)
saveRDS(out,'s1.rds')

## Scenario 2: Predominant Community (60%)
## Using tpar = c(0.21,0.83,0.91,0.67,1.08,1.2,1.2,1)
tpar <- matrix(rep(c(0.21,0.83,0.91,0.67,1.08,1.2,1.2,1),25),nrow = 25, byrow = T)
out <- parApply(cl,tpar,1,datagen,hhno1=500)
saveRDS(out,'s2.rds')

## Scenario 3: Predominant Household (60%)
## Using tpar = c(0.8,0.3,0.38,0.2,0.5,1.2,1.2,1)
tpar <- matrix(rep(c(0.8,0.3,0.38,0.2,0.5,1.2,1.2,1),25),nrow = 25, byrow = T)
out <- parApply(cl,tpar,1,datagen,hhno1=1000)
saveRDS(out,'s3.rds')


## Scenario 4: Extremely Predominant Household (Over 80%)
## Using tpar = c(1.2,0.1,0.15,0.05,0.2,1.2,1.2,1)
tpar <- matrix(rep(c(1.2,0.1,0.15,0.05,0.2,1.2,1.2,1),25),nrow = 25, byrow = T)
out <- parApply(cl,tpar,1,datagen,hhno1=2000)
saveRDS(out,'s4.rds')
stopCluster(cl)

