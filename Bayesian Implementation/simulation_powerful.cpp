
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
double like2(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector h1, 
             NumericVector h2, NumericVector h3, NumericVector h4, NumericVector h5,NumericVector c1,
             NumericVector c2, NumericVector c3, int index, int init, IntegerVector edge, 
             NumericVector par,double rn){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x3 = data["infect"];
   IntegerVector x4 = data["nopar"];
   IntegerVector h1 = data["crowd"];
   IntegerVector h2 = data["sleep"];
   IntegerVector h3 = data["time"];
   IntegerVector h4 = data["severity"];
   IntegerVector h5 = data["lag"];
   IntegerVector c1 = data["burden"];
   Integervector c2 = data["ses"];
   IntegerVector c3 = data["cough"];
   */
  double adult = x1[index];
  double female = x2[index];
  double infect = x3[index];
  double nopar = x4[index];
  double crowd = h1[index];
  double sleep = h2[index];
  double time = h3[index];
  double burden = c1[index];
  double ses = c2[index];
  double cough = c3[index];
  double prob = 0;
  int k = 0;
  if (infect == 2){
    double lag = h5[index];
    for (int i=init;i<init+nopar-1;++i){
      double severity = h4[i];
      double lag1 = h5[i];
      double lambdah = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-crowd)+crowd*par[3])*
        ((1-sleep)+sleep*par[4])*((1-time)+time*par[5])*((1-severity)+severity*par[6]);
      double probh = log((1-exp(-lambdah))*R::dgamma(lag-lag1,1.46,1/0.657,FALSE))*edge[k]+(1-edge[k])*log(1-(1-exp(-lambdah))*R::dgamma(lag-lag1,1.46,1/0.657,FALSE));
      prob += probh;
      ++k;
    }
    double lambdac = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-burden)+burden*par[7])*
      ((1-ses)+ses*par[8])*((1-cough)+cough*par[9]);
    prob += log((1-exp(-lambdac))*R::dgamma(rn,1.46,1/0.657,FALSE))*edge[k]+(1-edge[k])*log(1-(1-exp(-lambdac))*R::dgamma(rn,1.46,1/0.657,FALSE));
  }
  if (infect == 1){
    for (int i=init;i<init+nopar-1;++i){
      double severity = h4[i];
      double lambdah = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-crowd)+crowd*par[3])*
        ((1-sleep)+sleep*par[4])*((1-time)+time*par[5])*((1-severity)+severity*par[6]);
      double probh = log((1-exp(-lambdah)))*edge[k]-(1-edge[k])*lambdah;
      prob += probh;
      ++k;
    }
    double lambdac = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-burden)+burden*par[7])*
      ((1-ses)+ses*par[8])*((1-cough)+cough*par[9]);
    prob += log((1-exp(-lambdac)))*edge[k]-(1-edge[k])*lambdac;
  }
  return prob;
}


// for case with nopar == 0 and infect == 0
// [[Rcpp::export]]
double like3(NumericVector x1, NumericVector x2, NumericVector x5, NumericVector h1,
             NumericVector h2, NumericVector h3, NumericVector h4, NumericVector c1, 
             NumericVector c2, NumericVector c3, int index, int init, NumericVector par){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x5 = data["nopar1"];Note: nopar1 is nonzero for non-TB cases in HHC data
   IntegerVector h1 = data["crowd"];
   IntegerVector h2 = data["sleep"];
   IntegerVector h3 = data["time"];
   IntegerVector h4 = data["severity"];
   IntegerVector c1 = data["burden"];
   Integervector c2 = data["ses"];
   IntegerVector c3 = data["cough"];
   */
  double adult = x1[index];
  double female = x2[index];
  double nopar1 = x5[index];
  double crowd = h1[index];
  double sleep = h2[index];
  double time = h3[index];
  double burden = c1[index];
  double ses = c2[index];
  double cough = c3[index];
  double prob;
  double lambdac = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-burden)+burden*par[7])*
    ((1-ses)+ses*par[8])*((1-cough)+cough*par[9]);
  prob = -lambdac;
  for (int i=init;i<init+nopar1-1;++i){
    double severity = h4[i];
    double lambdah = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-crowd)+crowd*par[3])*
      ((1-sleep)+sleep*par[4])*((1-time)+time*par[5])*((1-severity)+severity*par[6]);
    double probh = -lambdah;
    prob += probh;
  }
  return prob;
} 

// The log-likelihood function like1 for the whole data
// [[Rcpp::export]]
double like1c(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector c1, 
              NumericVector c2, NumericVector c3, int index, NumericVector par){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x3 = data["infect"];
   IntegerVector c1 = data["burden"];
   Integervector c2 = data["ses"];
   IntegerVector c3 = data["cough"];
   */
  double adult = x1[index];
  double female = x2[index];
  double infect = x3[index];
  double burden = c1[index];
  double ses = c2[index];
  double cough = c3[index];
  double prob;
  double lambdac = par[0]*((1-adult)+adult*par[1])*((1-female)+female*par[2])*((1-burden)+burden*par[7])*
    ((1-ses)+ses*par[8])*((1-cough)+cough*par[9]);;
  prob = log(1-exp(-lambdac))*(infect>0)-lambdac*(infect==0);
  return prob;
}

// The log-likelihood function like for the whole data
// [[Rcpp::export]]
double likec(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5,
             NumericVector h1,NumericVector h2, NumericVector h3, NumericVector h4, NumericVector h5,
             NumericVector c1,NumericVector c2, NumericVector c3, NumericVector d1, IntegerVector t,
             IntegerVector edge, NumericVector par,int size,double rn){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x3 = data["infect"];
   IntegerVector x4 = data["nopar"];
   IntegerVector x5 = data["nopar1"];
   IntegerVector h1 = data["crowd"];
   IntegerVector h2 = data["sleep"];
   IntegerVector h3 = data["time"];
   IntegerVector h4 = data["severity"];
   IntegerVector h5 = data["lag"];
   IntegerVector c1 = data["burden"];
   Integervector c2 = data["ses"];
   IntegerVector c3 = data["cough"];
   IntegerVector d1 = data["hhc"];
   */
  double result = 0;
  int noedge = 0;
  for (int i=0;i<size;++i){
    if (d1[i]==0){
      result += like1c(x1,x2,x3,c1,c2,c3,i,par);
    }else {
      if((x4[i]==0) & (x3[i]==2)){
        result += like1c(x1,x2,x3,c1,c2,c3,i,par);
      }else if(x4[i]>0){
        int start_ind = initial(i,t);
        IntegerVector s = edge[Range(noedge,noedge+x4[i]-1)];
        noedge += x4[i];
        result +=like2(x1,x2,x3,x4,h1,h2,h3,h4,h5,c1,c2,c3,i,start_ind,s,par,rn);
      }else if(x3[i]==0){
        int start_ind = initial(i,t);
        result +=like3(x1,x2,x5,h1,h2,h3,h4,c1,c2,c3,i,start_ind,par);
      }
    } 
  }
  return result;
}

// the like function, compute the likelihood for all individuals
// [[Rcpp::export]]
double like(NumericVector x1, NumericVector x2, NumericVector x3, NumericVector x4, NumericVector x5,
            NumericVector h1,NumericVector h2, NumericVector h3, NumericVector h4, NumericVector h5,
            NumericVector c1,NumericVector c2, NumericVector c3, IntegerVector t,IntegerVector edge, 
            NumericVector par,int size, double rn){
  /* IntegerVector x1 = data["adult"];
   IntegerVector x2 = data["female"];
   IntegerVector x3 = data["infect"];
   IntegerVector x4 = data["nopar"];
   IntegerVector x5 = data["nopar1"];
   IntegerVector h1 = data["crowd"];
   IntegerVector h2 = data["sleep"];
   IntegerVector h3 = data["time"];
   IntegerVector h4 = data["severity"];
   IntegerVector h5 = data["lag"];
   IntegerVector c1 = data["burden"];
   Integervector c2 = data["ses"];
   IntegerVector c3 = data["cough"];
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
      result += like2(x1,x2,x3,x4,h1,h2,h3,h4,h5,c1,c2,c3,i,start_ind,s,par,rn);
    }else if(x3[i]==0){
      int start_ind = initial(i,t);
      result += like3(x1,x2,x5,h1,h2,h3,h4,c1,c2,c3,i,start_ind,par);
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
  NumericVector x3 = data["infect"];
  NumericVector x4 = data["nopar"];
  NumericVector x5 = data["nopar1"];
  NumericVector h1 = data["crowd"];
  NumericVector h2 = data["sleep"];
  NumericVector h3 = data["time"];
  NumericVector h4 = data["severity"];
  NumericVector h5 = data["lag"];
  NumericVector c1 = data["burden"];
  NumericVector c2 = data["ses"];
  NumericVector c3 = data["cough"];
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
  NumericMatrix output (iter,10);
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
    for (int j=0;j<10;++j){
      if (j == 0){
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma1));
        if (newpar[j] < 1) {
        ratio = (newpar[j]/oldpar[j])*exp(like(x1,x2,x3,x4,x5,h1,h2,h3,h4,h5,c1,c2,c3,total,path,newpar,size,rn1)
                                            -like(x1,x2,x3,x4,x5,h1,h2,h3,h4,h5,c1,c2,c3,total,path,oldpar,size,rn1));
        } else {
        ratio = 0;
        }
        decision = {ratio, 1};
        output(i,j) = (R::runif(0,1)<min(decision))?newpar[j]:oldpar[j];
        oldpar[j] = output(i,j);
      } else {
        newpar = clone(oldpar);
        newpar[j] = exp(R::rnorm(log(oldpar[j]),sigma2));
        ratio = (R::dnorm(log(newpar[j]),0,1,false))/(R::dnorm(log(oldpar[j]),0,1,false))*
          exp(like(x1,x2,x3,x4,x5,h1,h2,h3,h4,h5,c1,c2,c3,total,path,newpar,size,rn1)-
          like(x1,x2,x3,x4,x5,h1,h2,h3,h4,h5,c1,c2,c3,total,path,oldpar,size,rn1));
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
        ratio = edge_prop(old_path,test_path)*exp(like2(x1,x2,x3,x4,h1,h2,h3,h4,h5,c1,c2,c3,k,start,test_path,para,rn1)
                  -like2(x1,x2,x3,x4,h1,h2,h3,h4,h5,c1,c2,c3,k,start,old_path,para,rn1));
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
