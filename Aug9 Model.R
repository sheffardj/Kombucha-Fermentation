library(deSolve)
library(naniar)
library(ggplot2)
library(numDeriv)
library(nanier)


#######################  Custom Optimize Function my.optim()  #######################
my.optim = function(par, fn, ..., loop=FALSE) {
  args = list(...)
  args$method = "BFGS"
  Last = NA
  if (is.null(args$control)) args$control = list()
  args$control$ndeps = rep(1e-6, length(par))
  hess = hessian(fn, par, method.args=list(eps=1e-7,d=1e-6,zero.tol=1e-12))
  e = eigen(hess)
  repeat {
    steps = 10/sqrt(abs(e$values))
    steps[is.infinite(steps)] = 1e-6
    args$fn = function(x) fn(par + c(e$vectors %*% x))
    args$par=rep(0, length(par))
    args$control$parscale = steps
    args$hessian=FALSE
    o = do.call(optim, args)
    o$par = par + c(e$vector %*% o$par)
    hess = hessian(fn, o$par, method.args=list(eps=1e-7,d=1e-6,zero.tol=1e-12))
    e = eigen(hess)
    o$hessian=hess
    if (all(e$values > 0) || !loop) break;
    if (is.na(Last)) cat(o$value, " ") else cat(Last - o$value, " ")
    if (is.na(Last) || Last > o$value + 1e-8) Last = o$value else break
    par = o$par
  }
  cat("value=", o$value, "\n")
  o
}
#######################  Model   #######################
Model3<-function(t,init,params){
  with(as.list(c(init,params)),{
    dS  = -a * S*Y/(1+alpha*S)
    dF  = a * S*Y /(1+alpha*S)   - b*Y * F/(1+beta*F) - c*B*F/(1+epsilon*F)
    dE  = 2*b*F*Y/(1+beta*F)   - d*E*B/(1+delta*E)
    dA  = d*E*B/(1+delta*E)
    #
    dGs = a * S*Y/(1+alpha*S)   - f*Gs*B/(1+gamma*Gs) - e*B*Gs*(1+lambda*Gs)
    dGa = f* Gs*B/(1+gamma*Gs)
    dC  = c*B*F/(1+epsilon*F)  + e*B*Gs/(1+lambda*Gs)
    #
    dY  = rho * (a*S*Y/(1+alpha*S) + b*Y*F/(1+beta*F))
    dB  = r * (c*B*F/(1+epsilon*F) + e*B*Gs*(1+lambda*Gs))
    list(c(dS,dF,dE,dA,dGs,dGa,dC,dY,dB))
  })
}

##############################################                 ##############################################  
##############################################  AGITATED CASE  ##############################################
##############################################                 ##############################################
  
#######################  Agitated Model Params   #######################
if(1){
  params.S = c(
  a = 0.1 ,             #S -> F + Gs
  b = 1 ,                #F -> E
  c = .1 ,              #F -> C
  d = 0.1 ,              #E -> A
  e = 0.001 ,           #Gs -> C
  f = 0.1 ,             #Gs -> Ga                                                                                         
  alpha   = 100 ,
  beta    = 100 ,
  gamma   = 100 ,
  delta   = 100 ,
  epsilon = 100 ,
  lambda  = 100 ,
  rho     = 1 ,
  r       = 1
)}
#######################  Init - Siever && Conv   #######################
conv2<-c(S=343.3,
         F=180.156,
         E=46.069,
         A=60.05196,
         Gs=180.15588,
         Ga=196.16, 
         C=162.1406, #g/mol "per glucose unit" - (?)
         Y=1,
         B=1
)
#Svr2<-matrix(NA,nrow=21,ncol=8)
#Svr2[,1:7]<-Svr[1:21,]
SvrC = c(0.2, 0, 0.2, 2.1, 0.2)
Svr2[1,8]<-SvrC[1]
for(i in c(10,20)){   Svr2[i+1,8]=SvrC[(i/10+1)]    } #stack in cellulose data to Sievers data frame
colnames(Svr2)<-c("t","S","F","E","A","Gs","Ga","C")
initS <- c(c(Svr2[1,-1]),Y=1,B=0.005)

#######################  logLik Sievers   #######################
logLik.S <- function(params){
  if(any(params < 0)) return(1e10)
  output.S<-ode(y = initS/conv2, times = Svr2[,1], func = Model3, parms = params.S, method = "ode45")  
  sol.S = output.S[,2:8]%*%diag(conv2[1:7])
  colnames(sol.S)=c("S","F","E","A","Gs","Ga","C")
  if(any(is.na(sol.S))) return(1e10)
  x=sum(c((sol.S[-1,]-Svr2[-1,-1])^2)/abs(sol.S[-1,]), na.rm = TRUE)
  if (any(c(sol.S)<0,na.rm=TRUE)) x*1000 else x
}
#######################  optim Sievers  #######################
for(i in 1:1){
  theta.S <- optim(par=params.S, fn=logLik.S, method="Nelder-Mead", 
                 control=list(maxit=200000, trace=1, parscale=pmax(0.001,params.S)))
  params.S= theta.S$par
}
    theta.S2 <- my.optim(par=params.S[1:14], fn=logLik.S,loop=2)
  params.S[1:14] = theta.S2$par        
#######################  SIEVER AGITATED MODEL PLOTS  #######################  
  if(1){
    output.S<-ode(y = initS/conv2, times = Svr2[,1], func = Model3, parms = params.S[1:14], method = "ode45")
    output.S[,2:7]<-output.S[,2:7]%*%diag(conv2[1:6])
    
    data2 = data.frame()
    
    for(i in 2:8){ data2 <- rbind(data2, data.frame(t=Svr2[,1],y=output.S[,i],var=colnames(output.S)[i],curve="fit"))}
    for(i in 2:7){ data2 <- rbind(data2, data.frame(t=Svr2[,1],y=Svr2[,i],var=colnames(Svr2)[i],curve="data"))}
    
    fig2 <- ggplot(data2, aes(x=t, y=y, col=curve,shape=curve, lty = curve)) + geom_line() +
      facet_wrap(~var, scale="free") +
      geom_point() + 
      scale_linetype_manual(values = c(fit=1,data=0)) +
      scale_shape_manual(values = c(fit=NA, data=16))
    print(fig2)
  }
  

##############################################                ##############################################  
##############################################  STATIC CASE   ##############################################
##############################################                ##############################################
  
#######################  Static Modle Params   #######################
  if(1){
    params.L = c(
      a = 0.5 ,             #S -> F + Gs
      b = 0.9 ,                #F -> E
      c = 2 ,              #F -> C
      d = 2 ,              #E -> A
      e = 0.5 ,           #Gs -> C
      f = 0.5 ,             #Gs -> Ga                                                                                          
      alpha   = 60 ,
      beta    = 10 ,
      gamma   = 1 ,
      delta   = 400 ,
      epsilon = 50 ,
      lambda  = 1000 ,
      rho     = 10 ,
      r       = 5,
      A0.L    = 0.1,       # 15
      Ga0.L   = 0.5,       # 16
      Y0.L    = 1,           # 17
      B0.L    = 1            # 18
    )}
#######################  Loncar Cellulose Data Setup  #######################
  if(0){
    L22_C = matrix(nrow=11, ncol=2,byrow=TRUE) 
    L22_C[,1] <- 0:10
    totalMass <- rowSums(Lon22[,-1])
    for(i in 1:11){
      L22_C[i,2] <- c(totalMass[1]-totalMass[i])
      colnames(L22_C) <- c("day","Cellulose")
    }
    Lon22 <- cbind(Lon22, C=L22_C[,2])
  }
#######################  Loncar init  #######################
 if(0){
   initL <- c(Lon22[1,2], #S
             Lon22[1,3], #F
             Lon22[1,4], #E
             A=params.L["A0.L"][[1]], #A to fit 
             Lon22[1,5], #Gs
             Ga=params.L["Ga0.L"][[1]], #Ga to fit
             C=0,
             Y=params.L["Y0.L"][[1]],  # to fit
             B=params.L["B0.L"][[1]])  # to fit
 }
#######################  logLik Loncar  #######################
logLik.L <- function(params) {
  initL <- c(Lon22[1,2], #S
             Lon22[1,3], #F
             Lon22[1,4], #E
             A=params.L["A0.L"][[1]], #A to fit 
             Lon22[1,5], #Gs
             Ga=params.L["Ga0.L"][[1]], #Ga to fit
             C=0,
             Y=params.L["Y0.L"][[1]],  # to fit
             B=params.L["B0.L"][[1]])  # to fit
  if(any(params < 0)) return(1e10)
  output.L<-ode(y = initL/conv2, times = Lon22[,1], func = Model3, parms = params.L[1:14], method = "ode45")
  sol.L = output.L[,2:8]%*%diag(conv2[1:7])
  colnames(sol.L)=c("S","F","E","A","Gs","Ga","C")
  if(any(is.na(sol.L))) return(1e10)
  TA22<-sol.L[,4]+sol.L[,6]
  s1=sum(c((sol.L[-1,c(-4,-6)]-Lon22[-1,c(-1,-6)])^2/abs(sol.L[-1,c(-4,-6)])), na.rm = TRUE)
  s2=sum(c((TA22[2:11]-Lon22[-1,6])^2/abs(TA22[2:11])), na.rm=TRUE)
  s =s1+s2
  if (any(c(sol.L)<0, na.rm=TRUE)) s*1000 else s
}
#######################  optim Loncar  #######################
  theta.L <- optim(par=params.L, fn=logLik.L, method="Nelder-Mead", 
                  control=list(maxit=200000, trace=1, parscale=pmax(0.001,params.L)))
  params.L= theta.L$par
  theta.L <- my.optim(par=params.L, fn=logLik.L)
  params.L= theta.L$par
#######################  LONCAR STATIC MODEL PLOT  #########
  output.L<-ode(y = initL/conv2, times = Lon22[,1], func = Model3, parms = params.L[1:14],method = "ode45")
  output.L[,2:8]<-output.L[,2:8]%*%diag(conv2[1:7])
  
  data3 = data.frame()
  
  for(i in c(2:4,6,8)){ data3 <- rbind(data3, data.frame(t=Lon22[,1],y=output.L[,i],var=colnames(output.L)[i],curve="fit"))}
  data3 <- rbind(data3, data.frame(t=Lon22[,1],y=output.L[,5]+output.L[,7],var="TA",curve="fit"))                    # Add TA to data3
  for(i in c(9,10)){ data3 <- rbind(data3, data.frame(t=Lon22[,1],y=output.L[,i],var=colnames(output.L)[i],curve="fit"))}
  for(i in c(2:5,7)){ data3 <- rbind(data3, data.frame(t=Lon22[,1],y=Lon22[,i],var=colnames(Lon22)[i],curve="data"))}
  data3 <- rbind(data3, data.frame(t=Lon22[,1],y=Lon22[,6],var=colnames(Lon22)[6],curve="data"))
  
  fig3 <- ggplot(data3, aes(x=t, y=y, col=curve,shape=curve, lty = curve)) + geom_line() +
    facet_wrap(~var, scale="free") +
    geom_point() + 
    scale_linetype_manual(values = c(fit=1,data=0)) +
    scale_shape_manual(values = c(fit=NA, data=16))
  print(fig3)


output.L