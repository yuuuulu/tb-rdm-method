
#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
IntegerVector edge_sampler(int x) {
  IntegerVector edge (x);
  IntegerVector ind = seq(0,x-1);
  IntegerVector ind1 = sample(ind,1);
  edge[ind1] = 1;
  return edge;
}


// [[Rcpp::export]]
IntegerVector edge_updater(IntegerVector edge) {
  IntegerVector edge1 = clone(edge);
  int no = sum(edge1);
  int size = edge1.size();
  IntegerVector ind = seq(0,size-1);
  if (no == 1){
    LogicalVector x = (edge1 == 0);
    IntegerVector ind0 = ind[x];
    IntegerVector select = sample(ind0,1);
    edge1[select] = 1;
  } else {
    IntegerVector select = sample(ind,1);
    int s = select[0];
    int k = (edge1[s]==1)?0:1;
    edge1[s] = k;
  }
  return edge1;
}



// [[Rcpp::export]]
double edge_prop(IntegerVector oldpath, IntegerVector newpath) {
  IntegerVector diff = newpath - oldpath;
  double size = oldpath.size();
  double b = sum(oldpath);
  int change = sum(diff);
  double r;
  if (change == 1){
   r = (size-b)/(b+1);
  } 
  if (change == -1){
   r = b/(size-b+1);
  }
  return r;
}



// [[Rcpp::export]]
double edge_weight(double number, double rm, IntegerVector oldpath, IntegerVector newpath) {
  IntegerVector diff = newpath - oldpath;
  int size = oldpath.size();
  double weight;
  IntegerVector ind = seq(0,size-1);
  int change = sum(diff);
  if (change == 1){
    LogicalVector compare = (diff == 1);
    IntegerVector loc = ind[compare];
    if (rm > 0){
    weight = (loc[0]==(size-1))?(1/number):1;
    } else {
    weight = (loc[0]==(size-1))?1:number;
    }
  } 
  if (change == -1){
    LogicalVector compare = (diff == -1);
    IntegerVector loc = ind[compare];
    if (rm > 0){
    weight = (loc[0]==(size-1))?number:1;
    } else {
    weight = (loc[0]==(size-1))?1:(1/number);
    }
  }
  return weight;
}



// x is a result of t=as.vector(table(dat$hhid))
// Find the starting index for the household of a case given his/her index
// [[Rcpp::export]]
int initial(int a,IntegerVector x){
  IntegerVector x1 = cumsum(x);
  IntegerVector k = ifelse(x1 < (a+1), 1, 0);
  int s = sum(k);
  int s1;
  if (s == 0){
    s1 = 0;
  } else {
    s1 = x1[s-1];
  }
  return s1;
}




// for cases with nopar > 0 
// The last edge is the community edge
// [[Rcpp::export]]
double like2(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector mun1, 
             NumericVector mun2, NumericVector mun3, NumericVector mun4, NumericVector h1,NumericVector i1,
             int index, int init, IntegerVector edge, NumericVector par,double rn){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x3 = data["infect"];
   IntegerVector mun1 = data["ca"];
   IntegerVector mun2 = data["se"];
   IntegerVector mun3 = data["vi"];
   IntegerVector mun4 = data["vv"];
   IntegerVector x4 = data["nopar"];
   IntegerVector h1 = data["crowd"];
   IntegerVector i2 = data["smear"];
   NumericVector i1 = data["lag"];
   IntegerVector i3 = data["culture"];
   */
  double adult = x1[index];
  double female = x2[index];
  double infect = x3[index];
  double ca = mun1[index];
  double se = mun2[index];
  double vi = mun3[index];
  double vv = mun4[index];
  double nopar = x4[index];
  double crowd = h1[index];
  double prob = 0;
  int k = 0;
  if (infect == 2){
    double lag = i1[index];
    for (int i=init;i<init+nopar-1;++i){
      double lag1 = i1[i];
      double lambdah = par[0]*((1-crowd)+crowd*par[5])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
      double probh = log((1-exp(-lambdah))*R::dgamma(lag-lag1,1.46,1/0.657,FALSE))*edge[k]+(1-edge[k])*log(1-(1-exp(-lambdah))*R::dgamma(lag-lag1,1.46,1/0.657,FALSE));
      prob += probh;
      ++k;
    }
    double lambdac = (ca*par[1]+se*par[2]+vi*par[3]+vv*par[4])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
    prob += log((1-exp(-lambdac))*R::dgamma(rn,1.46,1/0.657,FALSE))*edge[k]+(1-edge[k])*log(1-(1-exp(-lambdac))*R::dgamma(rn,1.46,1/0.657,FALSE));
  }
  if (infect == 1){
    for (int i=init;i<init+nopar-1;++i){
      double lambdah = par[0]*((1-crowd)+crowd*par[5])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
      double probh = log((1-exp(-lambdah)))*edge[k]-(1-edge[k])*lambdah;
      prob += probh;
      ++k;
    }
    double lambdac = (ca*par[1]+se*par[2]+vi*par[3]+vv*par[4])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
    prob += log((1-exp(-lambdac)))*edge[k]-(1-edge[k])*lambdac;
  }
  return prob;
}


// for case with nopar == 0 and infect == 0
// [[Rcpp::export]]
double like3(NumericVector x1, NumericVector x2, NumericVector x5, NumericVector mun1,
             NumericVector mun2, NumericVector mun3, NumericVector mun4, NumericVector h1,
             int index, int init, NumericVector par){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector mun1 = data["ca"];
   IntegerVector mun2 = data["se"];
   IntegerVector mun3 = data["vi"];
   IntegerVector mun4 = data["vv"];
   IntegerVector x5 = data["nopar1"];Note: nopar1 is nonzero for non-TB cases in HHC data
   IntegerVector h1 = data["crowd"];
   IntegerVector i2 = data["smear"];
   IntegerVector i3 = data["culture"];
   */
  double adult = x1[index];
  double female = x2[index];
  double ca = mun1[index];
  double se = mun2[index];
  double vi = mun3[index];
  double vv = mun4[index];
  double nopar1 = x5[index];
  double crowd = h1[index];
  double prob;
  double lambdac = (ca*par[1]+se*par[2]+vi*par[3]+vv*par[4])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
  prob = -lambdac;
  for (int i=init;i<init+nopar1-1;++i){
    double lambdah = par[0]*((1-crowd)+crowd*par[5])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
    double probh = -lambdah;
    prob += probh;
  }
  return prob;
} 

// The log-likelihood function like1 for the whole data
// [[Rcpp::export]]
double like1c(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector mun1,
              NumericVector mun2, NumericVector mun3, NumericVector mun4, int index, NumericVector par){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x3 = data["infect"];
   IntegerVector mun1 = data["ca"];
   IntegerVector mun2 = data["se"];
   IntegerVector mun3 = data["vi"];
   IntegerVector mun4 = data["vv"];
   */
  double adult = x1[index];
  double female = x2[index];
  double ca = mun1[index];
  double se = mun2[index];
  double vi = mun3[index];
  double vv = mun4[index];
  double infect = x3[index];
  double prob;
  double lambdac = (ca*par[1]+se*par[2]+vi*par[3]+vv*par[4])*((1-adult)+adult*par[6])*((1-female)+female*par[7]);
  prob = log(1-exp(-lambdac))*(infect>0)-lambdac*(infect==0);
  return prob;
}

// The log-likelihood function like for the whole data
// [[Rcpp::export]]
double likec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5,
             NumericVector mun1,NumericVector mun2, NumericVector mun3, NumericVector mun4, NumericVector h1,
             NumericVector i1, NumericVector d1, IntegerVector t,IntegerVector edge, NumericVector par,int size,double rn){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector mun1 = data["ca"];
   IntegerVector mun2 = data["se"];
   IntegerVector mun3 = data["vi"];
   IntegerVector mun4 = data["vv"];
   IntegerVector x3 = data["infect"];
   IntegerVector x4 = data["nopar"];
   IntegerVector x5 = data["nopar1"];
   IntegerVector h1 = data["crowd"];
   NumericVector i1 = data["lag"];
   IntegerVector d1 = data["hhc"];
   */
  double result = 0;
  int noedge = 0;
  for (int i=0;i<size;++i){
    if (d1[i]==0){
      result += like1c(x1,x2,x3,mun1,mun2,mun3, mun4, i, par);
    }else {
      if((x4[i]==0) & (x3[i]==2)){
        result += like1c(x1,x2,x3,mun1,mun2,mun3, mun4, i, par);
      }else if(x4[i]>0){
        int start_ind = initial(i,t);
        IntegerVector s = edge[Range(noedge,noedge+x4[i]-1)];
        noedge += x4[i];
        result +=like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, i, start_ind, s, par,rn);
      }else if(x3[i]==0){
        int start_ind = initial(i,t);
        result +=like3(x1,x2,x5,mun1,mun2,mun3,mun4,h1, i, start_ind, par);
      }
    } 
  }
  return result;
}

// the like function, compute the likelihood for all individuals
// [[Rcpp::export]]
double like(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5,
            NumericVector mun1,NumericVector mun2, NumericVector mun3, NumericVector mun4, NumericVector h1,
            NumericVector i1, IntegerVector t,IntegerVector edge, NumericVector par,int size, double rn){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector mun1 = data["ca"];
   IntegerVector mun2 = data["se"];
   IntegerVector mun3 = data["vi"];
   IntegerVector mun4 = data["vv"];
   IntegerVector x3 = data["infect"];
   IntegerVector x4 = data["nopar"];
   IntegerVector x5 = data["nopar1"];
   IntegerVector h1 = data["crowd"];
   NumericVector i1 = data["lag"];
   IntegerVector i2 = data["smear"];
   IntegerVector i3 = data["culture"];
   */
  double result = 0;
  int noedge = 0;
  for (int i=0;i<size;++i){
    if ((x4[i]==0) & (x3[i]==2)){
      continue;
    }else if(x4[i]>0){
      int start_ind = initial(i,t);
      IntegerVector s = edge[Range(noedge,noedge+x4[i]-1)];
      noedge += x4[i];
      result +=like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, i, start_ind, s, par,rn);
    }else if(x3[i]==0){
      int start_ind = initial(i,t);
      result +=like3(x1,x2,x5,mun1,mun2,mun3,mun4,h1, i, start_ind, par);
    }
  }
  return result;
}



// Gather all pieces and create the MCMC function for TB household contact study
// Don't have thinning and burn-in parameters set up, but can do that in post-processing
// MCMC without weights on updating the transmission paths
// [[Rcpp::export]]
NumericMatrix mcmc(DataFrame data, int iter, NumericVector par0, double sigma1, double sigma2){
  NumericVector x1 = data["adult"];
  NumericVector x2 = data["female"];
  NumericVector mun1 = data["ca"];
  NumericVector mun2 = data["se"];
  NumericVector mun3 = data["vi"];
  NumericVector mun4 = data["vv"];
  NumericVector x3 = data["infect"];
  NumericVector x4 = data["nopar"];
  NumericVector x5 = data["nopar1"];
  NumericVector h1 = data["crowd"];
  NumericVector i1 = data["lag"];
  IntegerVector hhid = data["hhid"];
  int size = data.nrows();
  int parno = sum(x4);
  IntegerVector total = table(hhid);
  IntegerVector path (parno);
  int noedge = 0;
  for (int i=0;i<size;++i){
    if (x4[i]>0){
      path[Range(noedge,noedge+x4[i]-1)] = edge_sampler(x4[i]);
      noedge += x4[i];
    } 
  }
  // Create output matrix
  NumericMatrix output (iter,8);
  output(0,_) = par0;
  NumericVector oldpar;
  NumericVector newpar;
  NumericVector decision;
  IntegerVector old_path;
  IntegerVector test_path;
  int start;
  double ratio;
  double rn1;
  // M-H sampler 
  for (int i=1;i<iter;++i){
    rn1 = R::rgamma(1.46,1/0.657);
    oldpar = output(i-1,_);
    for (int j=0;j<8;++j){
      if (j<5){
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma1));
        ratio = (newpar[j]/oldpar[j])*exp(like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,newpar,size,rn1)
                                            -like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,oldpar,size,rn1));
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      } else {
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma2));
        ratio = (R::dnorm(log(newpar[j]),0,1,false))/(R::dnorm(log(oldpar[j]),0,1,false))*
          exp(like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,newpar,size,rn1)-
          like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,oldpar,size,rn1));
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      }
    }
    NumericVector para = output(i,_);
    // Update the digraph
    noedge = 0;
    for (int k=0;k<size;++k){
      if (x4[k]>0){
        start = initial(k,total);
        old_path = path[Range(noedge,noedge+x4[k]-1)];
        test_path = edge_updater(old_path);
        ratio = edge_prop(old_path,test_path)*exp(like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, test_path, para,rn1)
                  -like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, old_path, para,rn1));
        decision = {ratio, 1};
        if (R::runif(0,1)<min(decision)) {path[Range(noedge,noedge+x4[k]-1)] = test_path;} 
        noedge += x4[k];
      } 
    }
  }
  
  return output;
  
}

// Weighted MCMC approach
// [[Rcpp::export]]
NumericMatrix wmcmc(DataFrame data, int iter, NumericVector par0, double sigma1, double sigma2, double rmean, double rstd){
  NumericVector x1 = data["adult"];
  NumericVector x2 = data["female"];
  NumericVector mun1 = data["ca"];
  NumericVector mun2 = data["se"];
  NumericVector mun3 = data["vi"];
  NumericVector mun4 = data["vv"];
  NumericVector x3 = data["infect"];
  NumericVector x4 = data["nopar"];
  NumericVector x5 = data["nopar1"];
  NumericVector h1 = data["crowd"];
  NumericVector i1 = data["lag"];
  IntegerVector hhid = data["hhid"];
  int size = data.nrows();
  int parno = sum(x4);
  IntegerVector total = table(hhid);
  IntegerVector path (parno);
  int noedge = 0;
  for (int i=0;i<size;++i){
    if (x4[i]>0){
      path[Range(noedge,noedge+x4[i]-1)] = edge_sampler(x4[i]);
      noedge += x4[i];
    } 
  }
  // Create output matrix
  NumericMatrix output (iter,8);
  output(0,_) = par0;
  NumericVector oldpar;
  NumericVector newpar;
  NumericVector decision;
  IntegerVector old_path;
  IntegerVector test_path;
  int start;
  double ratio;
  double rn1;
  double rn2;
  // M-H sampler 
  for (int i=1;i<iter;++i){
    rn1 = R::rgamma(1.46,1/0.657);
    rn2 = exp(R::rnorm(rmean,rstd));
    oldpar = output(i-1,_);
    for (int j=0;j<8;++j){
      if (j<5){
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma1));
        ratio = (newpar[j]/oldpar[j])*exp(like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,newpar,size,rn1)
                                            -like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,oldpar,size,rn1));
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      } else {
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma2));
        ratio = (R::dnorm(log(newpar[j]),0,1,false))/(R::dnorm(log(oldpar[j]),0,1,false))*
          exp(like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,newpar,size,rn1)-
          like(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,total,path,oldpar,size,rn1));
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      }
    }
    NumericVector para = output(i,_);
    // Update the digraph
    noedge = 0;
    for (int k=0;k<size;++k){
      if (x4[k]>0){
        start = initial(k,total);
        old_path = path[Range(noedge,noedge+x4[k]-1)];
        test_path = edge_updater(old_path);
        if (x3[k]==1) {
          ratio = edge_weight(rn2,rmean,old_path,test_path)*edge_prop(old_path,test_path)*
            exp(like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, test_path, para,rn1)
                  -like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, old_path, para,rn1));
        } 
        if (x3[k]==2) {
          ratio = edge_prop(old_path,test_path)*exp(like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, test_path, para,rn1)
                        -like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, old_path, para,rn1));
        }
        decision = {ratio, 1};
        if (R::runif(0,1)<min(decision)) {path[Range(noedge,noedge+x4[k]-1)] = test_path;} 
        noedge += x4[k];
      } 
    }
  }
  
  return output;
  
}


// Make a MCMC function for the whole data, not just the household contact study
// Don't have thinning and burn-in parameters set up, but can do that in post-processing
// [[Rcpp::export]]
NumericMatrix mcmc1(DataFrame data, int iter, NumericVector par0, double sigma1, double sigma2){
  NumericVector x1 = data["adult"];
  NumericVector x2 = data["female"];
  NumericVector mun1 = data["ca"];
  NumericVector mun2 = data["se"];
  NumericVector mun3 = data["vi"];
  NumericVector mun4 = data["vv"];
  NumericVector x3 = data["infect"];
  NumericVector x4 = data["nopar"];
  NumericVector x5 = data["nopar1"];
  NumericVector h1 = data["crowd"];
  NumericVector d1 = data["hhc"];
  NumericVector i1 = data["lag"];
  IntegerVector hhid = data["hhid"];
  int size = data.nrows();
  int parno = sum(x4);
  IntegerVector total = table(hhid);
  IntegerVector path (parno);
  int noedge = 0;
  for (int i=0;i<size;++i){
    if (x4[i]>0){
      path[Range(noedge,noedge+x4[i]-1)] = edge_sampler(x4[i]);
      noedge += x4[i];
    } 
  }
  // Create output matrix
  NumericMatrix output (iter,8);
  output(0,_) = par0;
  NumericVector oldpar;
  NumericVector newpar;
  NumericVector decision;
  IntegerVector old_path;
  IntegerVector test_path;
  int start;
  double ratio;
  double rn1;
  // M-H sampler 
  for (int i=1;i<iter;++i){
    rn1 = R::rgamma(1.46,1/0.657);
    oldpar = output(i-1,_);
    for (int j=0;j<8;++j){
      if (j<5){
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma1));
        ratio = (newpar[j]/oldpar[j])*exp(likec(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,d1,total,path,newpar,size,rn1)
                                            -likec(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,d1,total,path,oldpar,size,rn1));
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      } else {
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma2));
        ratio = (R::dnorm(log(newpar[j]),0,1,false))/(R::dnorm(log(oldpar[j]),0,1,false))*
          exp(likec(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,d1,total,path,newpar,size,rn1)-
          likec(x1,x2,x3,x4,x5,mun1,mun2,mun3,mun4,h1,i1,d1,total,path,oldpar,size,rn1));
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      }
    }
    NumericVector para = output(i,_);
    // Update the digraph
    noedge = 0;
    for (int k=0;k<size;++k){
      if (x4[k]>0){
        start = initial(k,total);
        old_path = path[Range(noedge,noedge+x4[k]-1)];
        test_path = edge_updater(old_path);
        ratio = edge_prop(old_path,test_path)*exp(like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, test_path, para,rn1)
                      -like2(x1,x2,x3,x4,mun1,mun2,mun3,mun4,h1,i1, k, start, old_path, para,rn1));
        decision = {ratio, 1};
        if (R::runif(0,1)<min(decision)) {path[Range(noedge,noedge+x4[k]-1)] = test_path;} 
        noedge += x4[k];
      } 
    }
  }
  
  return output;
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R


*/
