#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
int findvalue(double k,NumericVector x){
  NumericVector k1 = {k};
  IntegerVector s = match(k1,x);
  int s1 = as<int>(s);
  return s1-1;
}

// [[Rcpp::export]]
int findcolumn(NumericMatrix m,double k){
  NumericVector m1 = as<NumericVector>(m);
  int s = findvalue(k,m1);
  int n = m.cols();
  int c = s/n;
  return c;
}

// [[Rcpp::export]]
NumericVector vrbinom(NumericVector prob){
  int n = prob.size();
  NumericVector output (n);
  for (int i=0;i<n;i++){
    output[i] = R::rbinom(1,prob[i]);
  }
  return output;
}


// [[Rcpp::export]]
IntegerVector dropcase(IntegerVector x, int c){
  LogicalVector select (x!=c);
  IntegerVector x1 = x[select];
  return x1;
}



// [[Rcpp::export]]
void takerow (NumericMatrix m,int row, IntegerVector col, NumericVector y){
  int n = col.size();
  for (int i=0;i<n;i++){
    m(row,col[i]) = y[i];
  }
}


// [[Rcpp::export]]
void takediag (NumericMatrix m, NumericVector y){
  int n = m.cols();
  for (int i=0;i<n;i++){
    m(i,i) = y[i];
  }
}


// [[Rcpp::export]]
NumericVector takecolumns (NumericMatrix m, IntegerVector col){
  int n = col.size();
  int r = m.rows();
  NumericMatrix m1(r,n);
  for (int i=0;i<n;i++){
    m1(_,i) = m(_,col[i]);
  }
  NumericVector m2 = as<NumericVector>(m1);
  return m2;
}

// [[Rcpp::export]]
IntegerVector append(IntegerVector x, int y){
  x.push_back(y);
  IntegerVector x1 = as<IntegerVector>(x);
  return x1;
}

// [[Rcpp::export]]
NumericVector hhcheck(NumericVector crowd, NumericVector female, NumericVector adult, NumericVector ca, NumericVector se,
                      NumericVector vi, NumericVector vv, int size, NumericVector tpar){
  double ca1 = as<double>(unique(ca));
  double se1 = as<double>(unique(se));
  double vi1 = as<double>(unique(vi));
  double vv1 = as<double>(unique(vv));
  double crowd1 = as<double>(unique(crowd));
  NumericMatrix transmat (size, size);
  NumericMatrix datemat (size, size);
  IntegerVector id = seq(0,size-1);
  NumericVector lambdac = (ca1*tpar[1]+se1*tpar[2]+vi1*tpar[3]+vv1*tpar[4])*((1-adult)+adult*tpar[6])*((1-female)+female*tpar[7]);
  NumericVector probc = 1-exp(-lambdac);
  NumericVector ltbi;
  ltbi = vrbinom(probc);
  NumericVector actb ;
  actb = rep(0,size);
  LogicalVector infected;
  infected = (ltbi>0);
  NumericVector proba;
  proba = rep(0.1,sum(infected));
  actb[infected] = vrbinom(proba);
  NumericVector trans;
  trans = ltbi + actb;
  trans[0] = 2;
  NumericVector date;
  date = rep(0,size);
  LogicalVector progressed;
  progressed = (trans==2);
  date[progressed] = Rcpp::runif(sum(progressed),0,10);
  takediag(transmat,trans);
  takediag(datemat,date);
  NumericVector date_nonzero;
  date_nonzero = date[progressed];
  double t = min(date_nonzero);
  LogicalVector active;
  active = (date==t);
  IntegerVector actind;
  actind = id[active];
  int sender = actind[0];
  IntegerVector receiver;
  receiver =  dropcase(id,sender);
  NumericVector adult1;
  NumericVector female1;
  NumericVector lambdah;
  NumericVector prob1;
  NumericVector result;
  LogicalVector stop1;
  LogicalVector stop2;
  NumericVector date1;
  LogicalVector selectdate;
  IntegerVector actind1;
  IntegerVector receiver1;
  while (actind.size()<=size){
    adult1 = adult[receiver];
    female1 = female[receiver];
    lambdah = tpar[0]*(crowd1*tpar[5]+(1-crowd1))*(adult1*tpar[6]+(1-adult1))*(female1*tpar[7]+(1-female1));
    prob1 = 1-exp(-lambdah);
    ltbi = vrbinom(prob1);
    infected = (ltbi>0);
    if(is_true(any(infected))){
      actb = rep(0,receiver.size());
      proba = rep(0.1,sum(infected));
      actb[infected] = vrbinom(proba);
      trans = ltbi + actb;
    } else {
      trans = ltbi; 
    }
    date = rep(0,receiver.size());
    progressed = (trans==2);
    if(is_true(any(progressed))){
      result = (t+Rcpp::rgamma(sum(progressed),1.46,1/0.657));
      date[progressed] = result;
    }
    takerow(transmat,sender,receiver,trans);
    takerow(datemat,sender,receiver,date);
    stop1 = (trans<2);
    date1 = takecolumns(datemat,receiver);
    stop2 = (date1==0);
    if(is_true(all(stop1))&is_true(all(stop2))) {break;}
    selectdate = (date1>0);
    date_nonzero = date1[selectdate];
    t = min(date_nonzero);
    sender = findcolumn(datemat,t);
    actind1 = clone(actind);
    actind = append(actind1,sender);
    if (actind.size() == size){break;}
    receiver1 = clone(receiver);
    receiver = dropcase(receiver1,sender);
  } 
  NumericVector output (size);
  for (int i = 0; i<size; i++){
    output[i] = max(transmat(_,i));
  }
  NumericVector output1 = output[Range(1,size-1)];
  return output1;
}
  
  
  
  
// [[Rcpp::export]]
double adcheck(DataFrame data, int iter, NumericVector tv, NumericVector tpar){
    NumericVector x1 = data["adult"];
    NumericVector x2 = data["female"];
    NumericVector mun1 = data["ca"];
    NumericVector mun2 = data["se"];
    NumericVector mun3 = data["vi"];
    NumericVector mun4 = data["vv"];
    NumericVector h1 = data["crowd"];
    IntegerVector hhid = data["hhid"];
    int size = tv.length();
    int nohh = max(hhid);
    IntegerVector total = table(hhid);
    IntegerVector total1 = total - 1;
    NumericMatrix output (iter,size);
    NumericVector value;
    NumericVector yes (size);
    NumericVector prob;
    LogicalVector infected;
    LogicalVector selecthh;
    NumericVector crowd;
    NumericVector adult;
    NumericVector female;
    NumericVector ca;
    NumericVector se;
    NumericVector vi;
    NumericVector vv;
    int iterator;
    NumericVector out (size);
    for (int i=0; i<iter; i++){
      iterator = 0;
      for (int j=0; j<nohh; j++){
        selecthh = (hhid==(j+1));
        crowd = h1[selecthh];
        adult = x1[selecthh];
        female = x2[selecthh];
        ca = mun1[selecthh];
        se = mun2[selecthh];
        vi = mun3[selecthh];
        vv = mun4[selecthh];
        out[Range(iterator,iterator+total1[j]-1)] = hhcheck(crowd,female,adult,ca,se,vi,vv,total[j],tpar);
        iterator += total1[j];
      }
      output(i,_) = clone(out);
    }
    for (int k = 0; k < size; k++){
      value = output(_,k);
      infected = (value>0);
      yes[k] = sum(infected);
    }
    prob = yes/iter;
    LogicalVector check1;
    LogicalVector check2;
    check1 = (prob==1);
    check2 = (prob==0);
    prob[check1] = 0.9999;
    prob[check2] = 0.0001;
    NumericVector loglike;
    loglike = tv*log(prob)+(1-tv)*log(1-prob);
    return sum(loglike);
  }





/*** R

*/
