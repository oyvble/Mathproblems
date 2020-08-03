#Binary triangle task: Maximize the number of zeros in the whole triangle

#What sequence of upper layer will maximize this?
#Function: (Binary input)^2 -> Binary output
#(0,0)=>1, (1,1)=>1, (0,1)=>0
#Upper layer gives values to the next layer in triangle

##Examples:
#Example with 2 layers:
#  [x11][x12]
#     [x21]
#
#x11=1,x12=0 gives x21=0
#There are 3 solutions giving the maximum of 2 zeroes

library(Rcpp)
sourceCpp(code='
 #include <Rcpp.h> // added this
 using namespace Rcpp;

  // [[Rcpp::export]]
 int countC(LogicalVector x) {
  int cc = 0;
  int nextN,i;

  while(true) {
    cc += x.length() - sum(x); //count number of FALSE
    nextN = x.length() - 1; //size of next layer
    if(nextN==0) return cc; //return when hit 0-layer

    LogicalVector nextx;
    for(i=0; i<nextN; i++) {
      nextx.push_back( x(i)==x(i+1) ); //new output 
    }
    x = nextx;
  }
 }
')

n = 3
init = sample(c(T,F),n,replace=T)
countC(init)


optimTree = function(n=4) {
  #n = 4 #width of problem
  ncomb = 2^n 
  #Function: (Binary input)^2 -> Binary output
  #(0,0)=>1, (1,1)=>1, (0,1)=>0
  out = c(TRUE,FALSE)
  outL = list()
  for(i in 1:n) outL[[i]] = out
  combs =  expand.grid(outL) #get outcomes
  nNum = rep(NA,ncomb)
  for(j in 1:ncomb) {
    newx = as.logical(combs[j,,drop=FALSE])
    nNum[j] = countC(newx)
  }
  maxNum = max(nNum)
  ind = maxNum==nNum
  return(list(maxinput=combs[ind,],maxval=maxNum))
}

optimTree(10)


#A method to perform MCMC
optimTreeMCMC = function(n=4,mcmciter=1000,mcmcreps=10,seed=1,verbose=F,deltab=0.001) {
  #n = 4 #width of problem
  set.seed(seed)
  treeSize = (n+1)*n/2
#  plot(0,0,ty="n",xlim=c(0,mcmciter),ylim=c(0,treeSize))
  
  maxInput = numeric()
  maxVal = 0
  
  for(rep in 1:mcmcreps) {
    nacc = 0 #number of accepts
    nNum = rep(NA,mcmciter)
    indFlip = sample(1:n,mcmciter,replace=T)
    u = runif(mcmciter)
    b = deltab #0.001 #/treeSize #temoerature
    x = sample(c(T,F),n,replace=T)
    nNum[1] = countC(x)
  
    #propose new solution:
    for(it in 2:mcmciter) {
      x2 = x
      x2[indFlip[it]] = !x2[indFlip[it]] #flip
      num2 = countC(x2) #get updated score
      Aprob = exp( b*(num2-nNum[it-1]) )
      if(  Aprob >  u[it]) { #accept
        nacc = nacc + 1
        x = x2 #accept proposal
        nNum[it] = num2 #update with new count
        
        if(num2>=maxVal) {
          maxVal =num2
          maxInput = rbind(x2)
        } else if(num2==maxVal) {
          maxInput = unique( rbind(maxInput,x2) )
        }
      } else {
        nNum[it] = nNum[it-1] #use prev
      }
      b = b + deltab #0.001 #increment
    }
    accr = nacc/mcmciter
    if(verbose) print(paste0("acc=",round(accr,3)," - max=",maxVal))   
    plot(1:mcmciter,nNum,ty="l")
  } #end for each iters
  return(list(maxinput=maxInput,maxval=maxVal))
}



#Hypothetical Formula for exact maximum:
getmaxn = function(maxn,verbose=FALSE) {
  dn = 2 #difference term: changed dyn. with recursive
  lastodd = FALSE #is difference even?
  maxvn = 2
  for(n in 2:(maxn-1)) {
    maxvn = maxvn + dn #add maximum
    
    iseven = dn%%2==0
    if(iseven || lastodd) dn = dn + 1 #add if dn was even or of last was also odd
    if(!iseven && !lastodd) lastodd = TRUE #set to odd
    if(iseven) lastodd = FALSE
    if(verbose) print(paste0("n=",n+1,": ",maxvn)) 
  }
  return(maxvn)
}

#50 OK
#51 NOT OK:MC=867 vs 884

n = 50 #insert number of layers in the triangle
trueMax = getmaxn(n,verbose=F);print(paste0("Suggestive max=",trueMax))
mcmc = optimTreeMCMC(n,mcmciter=1000,mcmcreps=100,seed=1,verbose=T,deltab=0.0004)

#Ultimate challenge:
#n=100: 3367 is max
#Best suggestion with MCMC: 2833

#$maxinput
#[,1] [,2]  [,3]  [,4]  [,5]  [,6] [,7]  [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19]
#x2 TRUE TRUE FALSE FALSE FALSE FALSE TRUE FALSE TRUE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE FALSE
#[,20] [,21] [,22] [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36] [,37]
#x2  TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
#[,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52] [,53] [,54] [,55]
#x2 FALSE FALSE  TRUE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE
#[,56] [,57] [,58] [,59] [,60] [,61] [,62] [,63] [,64] [,65] [,66] [,67] [,68] [,69] [,70] [,71] [,72] [,73]
#x2 FALSE  TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE
#[,74] [,75] [,76] [,77] [,78] [,79] [,80] [,81] [,82] [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90] [,91]
#x2 FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE
#[,92] [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100]
#x2 FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE   TRUE
#$maxval
#[1] 2728