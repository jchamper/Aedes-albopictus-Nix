#!/bin/RScript

rm(list=ls())
wd<-"/Users/yuehan/Documents/Nix model/Final Model/8.0" #set working directory
setwd(wd)

#Packages----
library(ggplot2)

#ML framework----
source("/Users/yuehan/Documents/Nix model/Final Model/8.0/Main.R")

#Global Variables----

CAGE1.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage1.txt" 
CAGE2.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage2.txt"
CAGE3.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage3.txt"
CAGE4.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage4.txt"

CAGE5.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage5.txt"
CAGE6.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage6.txt"
CAGE7.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage7.txt"

CAGE8.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage8.txt"
CAGE9.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage9.txt"
CAGE10.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage10.txt"
CAGE11.PATH<-"/Users/yuehan/Documents/Nix model/Final Model/Data/Cage11.txt"

DR.a.X<-0.44 #drive conversion rate
DRR.a.X<-0.0 #germline resistance formation rate

#ML setting
DIST.ab<-0.01 #convergence distance for CI determination, all cages 
A.ab<-0.05 #type 1 error, all cages  P<0.05
 
D.LIM<-9.99 #add to 1 for upper limit for drive fitness estimate in ML


skip_mark1=c() #transitions excluded in likelihood calculation, cage1

TOSTORE.ab<-"/Users/yuehan/Documents/Nix model/Final Model/Data/OUTPUT" #store ML output, all cages

#Additional Functions----
#1. read_emp_data 
#x = path to empirical cage data (.txt)
#return = data frame formatted for ML approach
read_emp_data<-function(x){
  y<-read.table(x)
  colnames(y)<-c("g","rm","wtm","f")   
  y$N<-NULL
  return(y)
}

#2. plot_emp_data
#x = output from read_emp_data()
#return = drive carrier frequency plot over time 
plot_emp_data<-function(x){
  y<-x
  y$rm<-x$rm/(x$rm+x$wtm+x$f)   #rm ratio
  y$wtm<-x$wtm/(x$rm+x$wtm+x$f) #wtm ratio
  y$f<-x$f/(x$rm+x$wtm+x$f) #f ratio
  y$g<-as.factor(y$g)   
  g<-ggplot(data=y,aes(x=g,y=r))+geom_line(col="red",group=1)+geom_point(col="red")
  g<-g+theme_minimal()+coord_cartesian(ylim = c(0,1),xlim=c(0,25))
  g<-g+scale_y_continuous("f(r)",breaks = c(0,0.25,0.5,0.75,1))
  g<-g+xlab("Generation")+theme(text=element_text(size = 15))
  g<-g+scale_x_discrete("Generation",breaks=seq(0,25,5))
  return(g)
}

#Load & EDA Data----

cage1<-read_emp_data(CAGE1.PATH)
cage2<-read_emp_data(CAGE2.PATH)
cage3<-read_emp_data(CAGE3.PATH)
cage4<-read_emp_data(CAGE4.PATH)

cage5<-read_emp_data(CAGE5.PATH)
cage6<-read_emp_data(CAGE6.PATH)
cage7<-read_emp_data(CAGE7.PATH)

cage8<-read_emp_data(CAGE8.PATH)
cage9<-read_emp_data(CAGE9.PATH)
cage10<-read_emp_data(CAGE10.PATH)
cage11<-read_emp_data(CAGE11.PATH)

plot_emp_data(cage1)
plot_emp_data(cage2)
plot_emp_data(cage3)
plot_emp_data(cage4)
plot_emp_data(cage5)
plot_emp_data(cage6)
plot_emp_data(cage7)
plot_emp_data(cage8)
plot_emp_data(cage9)
plot_emp_data(cage10)
plot_emp_data(cage11)

#ML cage All(joint analysis)
obs.data<-list()
obs.data[[1]]<-cage1
obs.data[[2]]<-cage2
obs.data[[3]]<-cage3
obs.data[[4]]<-cage4
obs.data[[5]]<-cage5[1:2,]
obs.data[[6]]<-cage6[1:2,]
obs.data[[7]]<-cage7[1:2,]
obs.data[[8]]<-cage8
obs.data[[9]]<-cage9
obs.data[[10]]<-cage10
obs.data[[11]]<-cage11

# obs.data<-list()
# obs.data[[1]]<-cage5[1:2,]
# obs.data[[2]]<-cage6[1:2,]
# obs.data[[3]]<-cage7[1:2,]
# obs.data[[4]]<-cage8
# obs.data[[5]]<-cage9
# obs.data[[6]]<-cage10
# obs.data[[7]]<-cage11

my.res.iter<-list() #set storage for individual runs

my.lim.iter <- set_limits(fd.lim = D.LIM) 

my.lim.iter[[1]][2] <- max(my.lim.iter[[1]][2], 0.1) 

# 打印检查
#print("Limits after adjustment:")
#print(my.lim.iter)
#print("Limits:")

#print(paste("Skip transition:",skip_mark1))

#ML estimate:
#print("START ML")
ml<-optim(my.lim.iter[[2]],logL,method = "L-BFGS-B",control = list(fnscale=-1),
          lower = my.lim.iter[[1]],upper = my.lim.iter[[3]],
          lN=obs.data,lmpath=M.path,
          lg=G,ldr = DR.a.X,ldrr=DRR.a.X,
          ls=skip_mark1)
#print(ml)
#print(paste("Maximum Log-Likelihood (lnL_max):", ml$value))

#AICc calculation:
#print("START AICc")
n.trans<-sum(unlist(lapply(obs.data,function(k) nrow(k)-1)))
#print(paste("Number of transistions:",n.trans))
aicc<-calculate_AICc(flnL = ml$value,p = length(my.lim.iter[[2]]),n = n.trans)
#print(paste("AICc:",aicc))

#CI estimate:
#print("START CI")
ci<-calculate_CI(fml = ml,fN = obs.data,fdist = DIST.ab,fa = A.ab,
                 fg = G,fmpath = M.path,
                 fdr = DR.a.X, fdrr = DRR.a.X,fs=skip_mark1)
#print("CI interval:")
print(ci[1])

my.res.iter[[1]]<-ml
my.res.iter[[2]]<-aicc
my.res.iter[[3]]<-ci
toStore<-paste(TOSTORE.ab,"ML-cages.rds",sep="")
saveRDS(my.res.iter,toStore)
rm(my.res.iter,ci,aicc,n.trans,ml,toStore,my.lim.iter)
print("------------------")

# --- Transition ---


all_cages <- list(
  cage1=cage1, cage2=cage2, cage3=cage3, cage4=cage4,
  cage5=cage5, cage6=cage6, cage7=cage7, cage8=cage8,
  cage9=cage9, cage10=cage10, cage11=cage11
)


transition_results <- list()

for (cage_name in names(all_cages)) {
  current_df <- all_cages[[cage_name]]
  n_rows <- nrow(current_df)
  
  # (t -> t+1)
  for (i in 1:(n_rows - 1)) {
    
    # transition
    
    obs_pair <- list(current_df[i:(i+1), ])
    
    trans_label <- paste0(cage_name, "_G", current_df$g[i], "-to-G", current_df$g[i+1])
    print(paste("Modeling Transition:", trans_label))
    
    
    my.lim.iter <- set_limits(fd.lim = D.LIM)
    my.lim.iter[[1]][2] <- max(my.lim.iter[[1]][2], 0.1) 
    
    # ML
    
    ml <- optim(my.lim.iter[[2]], logL, method = "L-BFGS-B",
                control = list(fnscale = -1),
                lower = my.lim.iter[[1]], upper = my.lim.iter[[3]],
                lN = obs_pair, lmpath = M.path,
                lg = G, ldr = DR.a.X, ldrr = DRR.a.X,
                ls = skip_mark1)
    
    # CI
    ci <- tryCatch({
      calculate_CI(fml = ml, fN = obs_pair, fdist = DIST.ab, fa = A.ab,
                   fg = G, fmpath = M.path,
                   fdr = DR.a.X, fdrr = DRR.a.X, fs = skip_mark1)
    }, error = function(e) return(NULL))
    
    # 存储
    transition_results[[trans_label]] <- list(
      cage = cage_name,
      gen_start = current_df$g[i],
      gen_end = current_df$g[i+1],
      ml = ml,
      ci = ci
    )
    #print("CI interval:")
    print(ci[1])
  }
}

saveRDS(transition_results, paste0(TOSTORE.ab, "All_Transitions_ML.rds"))