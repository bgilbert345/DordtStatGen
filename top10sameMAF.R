
library(Matrix)
source("continuous_simulator")


#define Qstat
top=function(x){tail(sort(x),causal)}

Qstat = function(gen, phen.perm){
  mu.hat = rep(mean(phen), length(phen))
  S = crossprod(gen, phen.perm-mu.hat)
  t=apply(abs(S),2,top)
  apply(t^2,2,sum)
}




non.settings=c(10,50,100,500,1000,5000,10000)
pow.trials = 10000
permutes = 1000
causal = 10 #number of causal variants
n = 1000 #number of individuals
beta = rep(0.75, causal) #effect size
maf = .01
s = 1 #model standard deviation
mu = 0 #model mean

pow.mat = matrix(nrow=length(non.settings),ncol=pow.trials)

for(k in 1:length(non.settings))
  
{
  
  noncausal= non.settings[k] #number of noncausal variants
  print(noncausal)
  ps = vector("numeric",length=pow.trials)
  
  for(j in 1:pow.trials){
    #simulate data
    sim = cont_sim(n, maf, beta, mu, s, causal,causal+noncausal)
    gen = Matrix(sim$gen,sparse=T)
    phen = sim$phen
    phen.perm = matrix(nrow=length(phen), ncol=permutes, replicate(permutes,sample(phen)))
    
        test.statistic = Qstat(gen, phen)
        #permutation test
        reps=Qstat(gen,phen.perm)
        ps[j]=mean(reps>=test.statistic)
    }
  pow.mat[k,]=ps
  
}
pow=vector("numeric") 
for (h in 1:length(pow.mat[,1])){
   pow=c(pow,mean(pow.mat[h,]<0.05))}


params=list(non.settings=non.settings, pow.trials=pow.trials, permutes=permutes, causal=causal,
n=n, beta=beta, maf=maf, s=s, mu=mu)

result=list(pow.mat=pow.mat, params=params)

save(result, file="top10sameMAF.Rdata")

EOF