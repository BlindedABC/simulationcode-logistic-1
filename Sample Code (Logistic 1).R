library(segmented)

#############   Functions: Proposed Method ########################
NR_logistic<-function(tau0)
{
  erro = FALSE;eta = 1; j = 1
  while (abs(eta) > 1e-11 & j<=500 & erro==FALSE)
  {
    tau0_1 = tau0
    #Step 1: fix the change-point and fit a logistic regression
    X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
    glm.inital <- glm(Y~X+ X_tau_1+Z,family="binomial",control = list(maxit = 25))
    theta0 = c(as.numeric(coef( glm.inital)),tau0_1)
    mu0 = predict(glm.inital, type="response")
    I_tau_1 = as.numeric(X > tau0_1)
    
    #Step 2: update the estimate of the change-point
    ### Calculate U and J ###
    Un =  -mean((Y-mu0)*I_tau_1)
    Jn =theta0[3]*mean(mu0*(1-mu0)*I_tau_1)
    tryCatch( {eta = Un/Jn},error=function(cond) {erro <<-  TRUE})
    if(erro ==TRUE| any(is.na(eta))) {eta=0}
    if (erro == FALSE)
    {tau_up = tau0 + eta # Update Tau
    tau0 = tau_up}
    j=j+1
  }
  tau.out = NA
  Test = (erro == FALSE && !is.na(eta) & max(tau0,na.rm=T)<max(X) & min(tau0,na.rm=T)>min(X))
  if ( Test )
  { tau.out = tau0}
  return(tau.out)
}

sd_logistic<- function(tau.out)
{
  erro=FALSE
  tau0_1 =  tau.out
  
  X_tau_1 =  as.numeric(X > tau0_1) * ( X - tau0_1 )
  I_tau_1 = as.numeric(X > tau0_1);
  
  glm.inital <- glm(Y~X+ X_tau_1+Z,family="binomial",
                    control = list(maxit = 25))
  theta0 = c(as.numeric(coef(glm.inital)),tau0_1)
  mu = predict(glm.inital, type="response")
  
  mu0 = mu*(1-mu)
  Jn1 =  c(mean(mu0),          mean(mu0*X),          mean(mu0*X_tau_1),         mean(mu0*Z),       
           mean(mu0*X),        mean(mu0*X^2),        mean(mu0*X*X_tau_1),       mean(mu0*X*Z), 
           mean(mu0*X_tau_1),  mean(mu0*X*X_tau_1),  mean(mu0*X_tau_1^2),       mean(mu0*X_tau_1*Z),
           mean(mu0*Z),        mean(mu0*X*Z),        mean(mu0*X_tau_1*Z),       mean(mu0*Z^2)
  )
  Jn1 = matrix( Jn1 , ncol=4,nrow=4)
  
  Jn2 = - theta0[3]*c(mean(mu0*I_tau_1),                                
                      mean(mu0*X*I_tau_1),                            
                      mean(mu0*X_tau_1),
                      mean(mu0*Z*I_tau_1)
  )
  Jn2 = matrix( Jn2 , ncol=1,nrow=4,byrow = T)
  Jn3 = t(Jn2)
  
  Jn4 = theta0[3]^2*mean(mu0*I_tau_1)
  Jn4 = matrix(Jn4,ncol=1,nrow=1)
  
  sd.out = NA
  tryCatch( { V = solve(Jn4 - Jn3%*%solve(Jn1)%*%Jn2)/(n)},
            error=function(cond) {erro <<- TRUE})
  
  if(erro==FALSE){  sd.out = sqrt(V)}
  return(sd.out)
}

############## Simulation Codes for Logistic 1 ##################
set.seed(123456)
#### true values beta0, beta1, beta2, tau, eta ####
theta = c(0.1, -5, 9, 0.4); eta = 0.2 
tau00= 0.4#runif(1,0.25,0.55) #initial values of the change-point
Nsim = 1000

n=200#sample size
i.seg= 1;i = 1;o= 0
tau_segest1=NULL;sd_segest1 =NULL;cov_seg1=0
tau_est1 = NULL;sd_est1 = NULL; cov1 = 0
####### Simulation Procedure ###########
while (o < Nsim)
{
  erro = FALSE
  ### generate simulation dataset ###
  Z = rnorm(n,0,1)
  a1=  rnorm(n,0,1);u = (a1+Z)/sqrt(2)
  X= pnorm(u,mean=0,sd=1) #X and Z are correlated and X~uniform(0,1)
  t=theta[1] + theta[2]*X + theta[3] * as.numeric(X > theta[4]) * ( X - theta[4] )+ eta*Z
  pr = 1/(1+exp(-t)) 
  Y = rbinom(n,1,pr)
  
  
  ### segmented method ###
  o1<-glm(Y~X+Z,family=binomial)
  tryCatch( {o.seg=segmented(o1,seg.Z =~X,psi=tau00,seg.control(toll = 1e-11, it.max = 500,n.boot=0,random=FALSE))},
            error=function(cond) {erro <<- TRUE})
  if (erro == FALSE &!is.null(o.seg$psi))
  {seg.out = o.seg$psi[2]
  seg.sd = o.seg$psi[3]
  tau_segest1[i.seg]=seg.out
  sd_segest1[i.seg]=seg.sd 
  cov_seg1 = cov_seg1  + as.numeric(seg.out-1.96*seg.sd <= theta[4] & seg.out+1.96*seg.sd >=theta[4])
  i.seg=i.seg+1}
  
  ### proposed method ###
  tau.out = NR_logistic(tau00)
  if (!is.na(tau.out))
  {  
    sd.out= sd_logistic(tau.out)
    if (!is.na( sd.out))
    {sd_est1[i] = sd.out
    tau_est1[i] = tau.out
    cov1 = cov1 + as.numeric(tau.out-1.96*sd.out <= theta[4] & tau.out+1.96*sd.out>=theta[4])
    i=i+1}}

  o=o+1
  if(o%%100==0){print(o)}
}

#Simulation Results of Logistic 1#
#segmented method
bias_tau1_seg = mean(tau_segest1) - theta[4] 
sd_tau1_seg = sd(tau_segest1)
mse_tau1_seg=sqrt(bias_tau1_seg^2+sd_tau1_seg ^2)
mean_sd1_seg = mean(sd_segest1)
c('Bias* MCSD* MSE* AVESE*',round(c(bias_tau1_seg,sd_tau1_seg,mse_tau1_seg,mean_sd1_seg )*1000,3),
  'CP%',cov_seg1/(i.seg-1)*100,#Convergence Probability
  'CR',(i.seg-1)/o*100#Convergence Rate of the Algorithm
)

#proposed method
bias_tau1 = mean(tau_est1) - theta[4] 
sd_tau1 = sd(tau_est1)
mse_tau1=sqrt(bias_tau1^2+sd_tau1 ^2)
mean_sd1 = mean(sd_est1)
c('Bias* MCSD* MSE* AVESE*',round(c(bias_tau1,sd_tau1,mse_tau1,mean_sd1 )*1000,3),
  'CP',cov1/(i-1)*100,#Convergence Probability
  'CR',(i-1)/o*100#Convergence Rate of the Algorithm
)


