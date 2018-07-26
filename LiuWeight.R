
library(Matrix)
source("continuous_simulator_Liu")
 
 
 #define Qstat
 Qstat = function(gen, phen.perm, p,mafs,a1,a2){
   mu.hat = rep(mean(phen), length(phen))
   w=dbeta(mafs,a1,a2)^2
   S = crossprod(gen, phen.perm-mu.hat)
   if (p==Inf){
     return(apply(w*abs(S),2,max))
   }
   Sp = w*(abs(S)^p)
   colSums(Sp)^(1/p)
 }
 
 non.settings=c(8,48,100,500,1000,5000,10000)
 pow.trials = 10000
 p=c(2,4,Inf)
 permutes = 1000
 causal = 8 #number of causal variants
 n = 2000 #number of individuals
 beta = rep(0.75, causal) #effect size
 maf1 = .001
 maf2 = .01
 s = 1 #model standard deviation
 mu = 0 #model mean
 a1=1 #beta distribution parameter1 for weighting
 a2=25 #beta distribution paramter2 for weighting
 
 pow.mat = matrix(nrow=length(non.settings), ncol=length(p))
 
 for(k in 1:length(non.settings))
   
 {
   
   noncausal= non.settings[k] #number of noncausal variants
   ps = matrix(ncol=length(p),nrow=pow.trials)
   
   for(j in 1:pow.trials){
     #simulate data
     sim = cont_sim(n, maf1,maf2, beta, mu, s, causal,causal+noncausal)
     gen = Matrix(sim$gen,sparse=T)
     phen = sim$phen
     phen.perm = matrix(nrow=length(phen), ncol=permutes, replicate(permutes,sample(phen)))
     mafs=colSums(gen)/n
     for(m in 1:length(p)){
       test.statistic = Qstat(gen, phen,p[m],mafs,a1,a2)
       #permutation test
       reps=Qstat(gen,phen.perm,p[m],mafs,a1,a2)
       ps[j,m]=mean(reps>=test.statistic)
     }
   }
   for(r in 1:length(p)){
     pow.mat[k,r] = mean(ps[,r]<.05)
   }
 }
 
 colnames(pow.mat) = sapply(p, toString)
 row.names(pow.mat) = sapply(non.settings, toString)
 
 params=list(non.settings=non.settings, pow.trials=pow.trials, p=p, permutes=permutes, causal=causal, n=n, beta=beta, maf1=maf1, maf2=maf2, s=s, mu=mu)
 
 result=list(pow.mat=pow.mat, params=params)
 
 save(result, file="LiuWeight.Rdata")
 
 
