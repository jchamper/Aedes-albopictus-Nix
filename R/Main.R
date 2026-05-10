#!/bin/RScript
#ML framework
# fixed ne

#Packages----
library(stringr)
library(ggplot2)

#Global Variables----
M.path<-"/Users/yuehan/Documents/Nix model/Final Model/7.0/Linked_frequency.csv"

#w=wildtype, drive locus
#d=drive, drive locus
#r=resistance, drive locus
G<-c("dwxx", "ddxx", "wwxy", # male genotypes, the dd genotypes won't actually exist, but we need them to calculate.
     "wwxx") # female genotypes

#Drive Functions----

#init_m
#x ... path to raw transition matrix (counts)
#y ... names of count columns
#return(z) ... expected offspring frequencies
init_m<-function(x,y){
  z<-read.csv(x,header = T)
  z<-cbind(z[,1:2],z[,c(3:6)]/rowSums(z[,c(3:6)]))
  colnames(z)<-c("m","f",y)
  return(z)
}

#init_genotype
#x ... names of genotypes
#return(y) ... empty named genotype vector
init_genotype<-function(x){
  y<-rep(0,6)
  names(y)<-x
  return(y)
}

#germ_drive
#v ... original offspring frequencies (length=12)
#x ... gene drive conversion rate
#y ... germ line resistance formation rate
germ_drive<-function(v,x,y){
  w<-v
  #print(w)
  w["dwxx"]<-v["dwxx"]*(1-x-y)
  
  #the frequency of wrxx, wrxy, wwxx, wwxy, rrxx and rrxy won't change here.
  #although dd genotypes don't naturally exist, we need to calculate how dw genotype works.
  
  w["ddxx"]<-v["ddxx"]+v["dwxx"]*x
  
  #v["ddxx"] should always be 0.
  
  return(w)
}

#mod_transition_m
#m ... transition matrix 

#g ... names of genotypes
#dr ... drive conversion rate
#drr ... drive resistance formation rate

#return(y) ... transition matrix with adapted offspring frequencies after germline + embryo cutting
mod_transition_m<-function(m,g,dr,drr){
  y<-m
  #Step 1: Germline cuts 
  for(i in c(1:nrow(y))){
    fy_iter<-y[i,] #store raw  transition probabilities + cross
    fmale_genotype<-init_genotype(x = g) #empty male genotype vector
    ffemale_genotype<-init_genotype(x = g) #empty female genotype vector
    
    fmale_genotype[fy_iter[["m"]]]<-1 #set male genotype
    ffemale_genotype[fy_iter[["f"]]]<-1 #set female genotype
    
    fmale_genotype<-germ_drive(fmale_genotype,x = dr ,y = drr) #drive conversion + resistance
    
    fmale_genotype<-fmale_genotype[1:3]
    # print("fmale_genotype")
    # print(fmale_genotype)
    ffemale_genotype<-ffemale_genotype[4]
    # print("ffemale_genotype")
    # print(ffemale_genotype)
    
    fnew_mates<-as.vector(outer(ffemale_genotype,fmale_genotype))  #all possible mating probs
    #print("fnew_mates")
    #print(fnew_mates)
    fnew_progeny<-fnew_mates*m[,(3:6)] #determine new weights
    y[i,(3:6)]<-colSums(fnew_progeny)/sum(fnew_progeny) #replace marginal probs for offspring
  }
  #print("mod_transition_m")
  #print(y)
  return(y)
}

#determine_fitness_coefficients
# fg ... names of genotypes
# fd ... drive fitness
#return(x) ... fitness vector for genotypes (n=12); rel weights
determine_fitness_coefficients<-function(fg,f){
  x<-rep(1,4)
  names(x)<-fg
  fd<-f["fd"]
  has_d<-grepl("d", names(x))
  
  for(i in 1:length(x)){
    current_cost<-1
    if(has_d[i]){
      # Individual has Drive allele: apply Drive fitness cost
      current_cost<-current_cost*fd
    }
    x[i]<-current_cost
  }
  return(x)
}

#geno_to_pheno
#x ... named vector of genotypes (m+f,length=8)
#return= vector of phenotype frequencies | x (rm,wtm,f)
geno_to_pheno<-function(x){
  fd_id<-grep("d",names(x)) #"drxx", "dwxx", "ddxx", # the dd genotypes won't actually exist, but we need them to calculate.
  rm<-sum(x[fd_id])
  wtm_id<-c("wwxy")
  wtm<-sum(x[wtm_id])
  f_id<-c("wwxx")
  f<-sum(x[f_id])
  y<-c(rm,wtm,f)
  y<-y/sum(y)
  names(y)<-c("rm","wtm","f")
  #print(y)
  
  return(y)
}

#pheno_to_geno
#g ... named vector of genotypes (m+f, length=8)
#p ... observed phenotypes (rm,wtm,f)
#return(fgcorr) ... vector of genotype frequencies, corrected by observed phenotypes (m+f,length=8)
pheno_to_geno<-function(g,p){
  fpobs<-p/sum(p)
  fpraw<-geno_to_pheno(g) #marginal sums of genotypes giving the same phenotype 
  fgcorr<-numeric(length = 4)
  names(fgcorr)<-names(g)
  fd_id<-grep("d",names(fgcorr))
  wtm_id<-c("wwxy")
  f_id<-c("wwxx")
  
  fgcorr[fd_id]<-g[fd_id]/fpraw["rm"]*fpobs["rm"]   #standardize 
  fgcorr[wtm_id]<-g[wtm_id]/fpraw["wtm"]*fpobs["wtm"]
  fgcorr[f_id]<-g[f_id]/fpraw["f"]*fpobs["f"]
  return(fgcorr)
}

#start_genotype .... function to start a genotype 
#x ... counts (rm,wtm,f)
#fg ... genotype names
#return(y)... genotype frequencies m+f (length=8)
start_genotype<-function(x,fg){
  x<-unlist(x)
  names(x)<-c("rm","wtm","f")
  #print(x)
  x<-x/sum(x)
  y<-numeric(length = 4)
  names(y)<-fg
  y["dwxx"]<-x["rm"] 
  y["wwxy"]<-x["wtm"]
  y["wwxx"]<-x["f"]
  #print(y)
  
  return(y)
}

#propagate_genotypes
#p ... parental genotypes (8)
#m ... transition matrix (modified)
#g ... names of genotypes
#f... fitness parameters (named!)
#return=expected genotype frequencies next generation (8)
propagate_genotypes<-function(p,m,g,f){
  eps <- 1e-12

  p[1:3] <- p[1:3] + eps
  p[4] <- p[4] + eps
  p <- p / sum(p)
  
  fmale<-p[1:3]/sum(p[1:3]) #standardize male freqs
  #print(fmale)
  ffemale<-p[4]/sum(p[4]) #standardize female freqs
  #print(ffemale)
  fmate<-determine_fitness_coefficients(fg = g,f= f)
  #print(fmate)
  
  fmale.coeff <- fmate[1:3]
  
  fmale.weighted <- fmale * fmale.coeff
  fmale.sum <- sum(fmale.weighted)
  
  
  fmale <- fmale.weighted / fmale.sum
  #print(fmale)
  fpairs<-as.vector(outer(ffemale,fmale))
  fmmate<-fpairs*m[,3:6]
  #print(m[,3:10])
  #print(fmmate)
  
  #fecundity selection
  ffec<-determine_fitness_coefficients(fg = g,f= f)
  ffec<-rep(ffec,length(fmale))
  fmfec<-ffec*fmmate
  #print(fmfec)
  
  #viability selection
  fvia<-1.0
  fmvia<-t(fvia*t(fmfec))
  #print(fmvia)
  
  #extend to males + females again 
  g<-colSums(fmvia)/sum(fmvia)
  g<-g/sum(g)
  return(g)
  
}

#ML Functions----

#calc_rho
#fp ... predicted frequencies (rm,wtm,f)
#fo ... observed frequencies (rm,wtm,f)
#fn ... effective population size
#return(r) ...  logL calculation
calc_rho <- function(fp, fo, fn) {
  
  eps <- 1e-10
  fp <- (fp + eps) / sum(fp + eps)
  
  
  fo.count <- fo * fn
  
  
  # r = Σ(n_i * log(p_i))
  r <- sum(fo.count * log(fp))
  
  return(r)
}

logL<-function(lf,lN,lmpath,lg,ldr,ldrr,ls){
  
  #neutral set up:
  
  l.ne<-1 #effective population size
  l.fd<-1
  
  l.ne<-lf["ne"] #always estimate Ne
  l.m<-init_m(lmpath,lg)
  l.mmod<-mod_transition_m(m = l.m,g = lg,dr = ldr,drr = ldrr) #modified transition matrix
  l.fd<-lf["fd"]
  
  l.params<-c(l.ne,l.fd)
  names(l.params)<-c("ne","fd")
  
  l<-0 #start logL calculation
  l.n.sets<-length(lN)
  for(j in c(1:l.n.sets)){
    l.Nj<-lN[[j]] #store counts
    l.Nfreq<-data.frame(rm=l.Nj$rm,wtm=l.Nj$wtm,f=l.Nj$f) #normalize counts
    l.Nfreq$n<-l.Nfreq$rm+l.Nfreq$wtm+l.Nfreq$f
    l.Nfreq$rm<-l.Nfreq$rm/l.Nfreq$n
    l.Nfreq$wtm<-l.Nfreq$wtm/l.Nfreq$n
    l.Nfreq$f<-l.Nfreq$f/l.Nfreq$n
    l.Nfreq$n<-NULL
    if (j==1){l.obs.geno<-start_genotype(x = l.Nfreq[1,],fg = lg)} #initialize
    
    for(i in c(2:nrow(l.Nj))){ 
      l.obs.pheno<-unlist(l.Nfreq[i,]) #observed phenotype
      
      l.pred.geno<-propagate_genotypes(p = l.obs.geno,m = l.mmod,g = lg,f = l.params) #expected genotype
      
      l.pred.pheno<-geno_to_pheno(l.pred.geno) #expected phenotype
      
      ifelse(i%in%ls, l<-l, l<-l+calc_rho(fp = l.pred.pheno,fo = l.obs.pheno,fn = l.ne))#logL one transition step
      l.obs.geno<-pheno_to_geno(p = l.obs.pheno,g = l.pred.geno) #offspring in j  = parents in j+1
    }
  }
  return(l)
}

set_limits <- function(fd.lim = 9.99) {
  
  y <- c(ne = 25, fd = 0.5)
  
  
  y.down <- c(ne = 25, fd = 0.5)
  
  
  y.up <- c(ne = 1000, fd = fd.lim)
  
  return(list(y.down, y, y.up))
}

#calculate_AICc
#flnL ... log likelihood
#p ... number of params
#n ... number of data points
#return(aicc)... corrected AIC for small sample size
calculate_AICc<-function(flnL,p,n){
  aic=2*p-2*flnL+c(2*p*p+2*p)/c(n-p-1)
  return(aic)
}

# calculate_CI
# fml ... ML result (optim(logL))
# fN  ... observed data
# fdist ... min. dist for CI search
# fa ... type 1 error

# fmax_iter
calculate_CI<-function(fml,fN,fdist,fa,fg,fmpath,fdr,fdrr,fs,fmax_iter=1000){

  fdist_i<-fdist
  fsig<-qchisq(1-fa,1) # calculate Chi^2 stat
  fLmax<-fml$value 
  y<-data.frame(par=fml$par)
  y$low_ci<-c(NA)
  y$up_ci<-c(NA)

  for(i in c(1:length(fml$par))){
    par_name <- names(fml$par)[i]

    fhard_max <- ifelse(i==1, 10000, 100) # ne_prop 
    fhard_min <- ifelse(i==1, 25, fdist) # ne_prop 

    fpara_dynamic<-fml$par # init
    fpara_test<-fpara_dynamic[i]

    
    fstep<-max(fdist, fpara_test * 0.01)

    current_iter<-0
    frun<-T

    while(frun){
      current_iter <- current_iter + 1

      # Step
      fpara_test_new <- max(fhard_min, fpara_test - fstep)

     
      if(fpara_test_new == fpara_test) {
        frun <- F
        y$low_ci[i]<-fpara_test
        
        break
      }

      fpara_dynamic[i]<-fpara_test_new 
      fl<-logL(lf = fpara_dynamic,lN = fN,lmpath = fmpath,lg = fg,ldr = fdr,ldrr = fdrr,ls=fs)
      fstat<-2*(fLmax-fl) # Chi^2

      if(fstat > fsig){ 
        if(2*fstep < fdist) { 
          frun<-F
          y$low_ci[i]<-fpara_test_new
          
        } else {
          fpara_test<-fpara_test 
          fstep<-fstep/2 
        }
      } else { 
        fpara_test<-fpara_test_new 
      }

      
      if(current_iter >= fmax_iter){
        frun<-F
        y$low_ci[i]<-fpara_test 
        
      }
    }

    fpara_dynamic<-fml$par # init upstream
    fpara_test<-fpara_dynamic[i]
    fstep<-max(fdist, fpara_test * 0.01) 

    current_iter<-0
    frun<-T

    while(frun){
      current_iter <- current_iter + 1

      
      fpara_test_new <- fpara_test + fstep

      
      if(fpara_test_new >= fhard_max) {
        frun <- F
        y$up_ci[i]<-fhard_max
        
        break
      }

      fpara_dynamic[i]<-fpara_test_new
      fl<-logL(lf = fpara_dynamic,lN = fN,lmpath = fmpath,lg = fg,ldr = fdr,ldrr = fdrr,ls=fs)
      fstat<-2*(fLmax-fl) # Chi^2

      if(fstat > fsig){ 
        if(2*fstep < fdist) { 
          frun<-F
          y$up_ci[i]<-fpara_test_new
          
        } else {
          fpara_test<-fpara_test 
          fstep<-fstep/2 
        }
      } else { 
        fpara_test<-fpara_test_new 
      }

      
      if(current_iter >= fmax_iter){
        frun<-F
        y$up_ci[i]<-fpara_test 
        
      }
    }

  }
  return(y)
}