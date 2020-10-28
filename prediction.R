library(rootSolve)

# This code is provided as part of the Supporting Information for the article
#     "Generic assembly patterns in large ecological communities" (2017)
#     M. Barbier, J.-F. Arnoldi, G. Bunin and M. Loreau.

# The function "prediction_from_matrix(r,A,D)" takes matrices for 
# growth rates r, interactions A and self-interactions D and returns
# a list containing the four reference parameters 
#    zeta
#    mu
#    sigma
#    gamma
# as well as all predicted quantities
#    biomass: total biomass of the system
#    sdBiomass: standard deviation of individual biomass
#    phi: fraction of survivors
#    survivors: number of survivors
#    simpsonD: Simpson diversity (1/Simpson index)
#    productivity: productivity ratio
#    variability: response to stochastic perturbation
#    press: response to random press perturbation (if infinite, signals multistability)
#
# OPTIONS
#    correlations
#       when set to TRUE (default), the function uses the extension of the reference model
#       that accounts for first-order correlations between the coefficients, and returns
#         cKa: measure of correlation between carrying capacities and interactions
#         sigma_row: measure of row-wise correlations in interactions
#    FR
#       If provided as a list (f,df) consisting of a function f(z) and its derivative f'(z)
#       computes the predictions for the model with functional response f(z)


measure_parameters <- function(r,A,D,PK=TRUE,correlations=TRUE,plotK=FALSE)
{
  K <-r/D
  alpha <- A / rep(D, each = nrow(A))
  S <- length(K)
  if( mean(alpha)>1){
    print("WARNING: Average interactions are strong (A/D>1). Predictions will not be accurate.")
  }
  else if (max(alpha)>1){
    if (max(alpha)>5){
      print("WARNING: Some very strong interactions (A/D>>1). Predictions may not be accurate.")
    }
    else {
      print("WARNING: Some strong interactions (A/D>1). Predictions may not be accurate.")
    }
  }
  avgK <- mean(K)
  stdK <- sd(K)
  histK <- hist(K,plot=FALSE)
  Kmin <- min(K)-5*stdK
  Kmax <- max(K)+5*stdK
  if(PK){
    PK <- function(x) {
      res <- exp(predict(smooth.spline(x=histK$mids, y=log(histK$density)),x)$y)
      return(res)
    }
  }
  else{
    PK <- function(x){dnorm(x,avgK,stdK)}
  }
  if(plotK){
    ks<-seq(Kmin-.5,Kmax+.5,0.01)
    plot(ks,PK(ks),type="l" )
    lines(ks,dnorm(ks,avgK,stdK),col="red")
    points(histK$mids,histK$density )
  }
  
  offdiag <- function(a){ a[row(a)!=col(a)] }
  offa <- offdiag(alpha)
  mu <-  (S-1)* mean( offa )
  sig <- sqrt(S-1)*sd(offa ) 
  gam <- cor(as.vector(offa),as.vector(offdiag(t(alpha))))
  
  if(correlations){
    corrKA <- (S-1)* mean(offdiag(rep(K, each = nrow(alpha))*(alpha)) ) -  avgK* mu
    sigrow <- sqrt(S-1)*sd(rowMeans(alpha) )
    sig <- sqrt(sig^2 - sigrow^2)
  }
  else
  { corrKA <-0
  sigrow <-0 }
  
  return( list(S=S,mu=mu,sig=sig,gam=gam, PK=list(avgK=avgK,stdK=stdK,dist=PK,Kmin=Kmin,Kmax=Kmax),correl=c(corrKA,sigrow) ))
}



prediction <- function(parameters,correlations=TRUE,FR=FALSE)
{
  
  PK <- parameters$PK
  if(!(identical(PK$dist,FALSE) & identical(FR,FALSE) )  ){
    #Make integral table for interpolation
    Kmin=PK$Kmin
    Kmax=PK$Kmax
    if(identical(FR,FALSE)){
      #Linear functional response by default
      FR <- list(f=function(x){x},df=function(x){1})
    }
    Kstep <- (Kmax-Kmin)/100
    KS <- seq(Kmin,Kmax,Kstep)
    Ktable <-  mapply(function(mom){
                  mapply(function(z0){
                    integrate( function(x){ FR$f(x)^mom * PK$dist(x+z0) },0, Inf )$value
                  }, KS ) }, c(0,1,2) )
    
    parameters$PK$table <-
      mapply(function(mom){function(K){ 
        K[K>Kmax]=Kmax
        K[K<Kmin]=Kmin
        approxfun(KS , Ktable[,mom] )(K) }
        }, c(1,2,3)
      )
  }

  z_parameters <- function(x,S,mu,sig,gam){
    meanz <- mu*x[1]*x[2]
    varz <- sig^2 * x[1] *x[3] 
    return(c(meanz,varz) )
  }
  gaussian_parameters <- function(x,S,mu,sig,gam,PK,correl){
    g <-gam*sig^2 * x[1]
    v <- x[4]
    u <- 1-g*v
    zparam <- z_parameters(x,S,mu,sig,gam)
    mean0 <- PK$avgK - zparam[1]
    corrKA <- correl[1]
    sigrow <- correl[2]
    var0 <- PK$stdK^2 +sigrow^2 * (S-1)*(x[1] *x[2])^2 - 2*x[1]*x[2]*corrKA + zparam[2]
    return(c(mean0/u,var0/u^2,u))
  }

  erfmom <- function(mom,mean0,var0){
    #If functional response is linear and all distributions are normal
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    if (var0<=0.001){var0<-0.001}
    xx<-mean0/sqrt(2*var0)
    mom0 <- .5 * (erf(xx)+1 )
    mom1 <- sqrt(var0/2/pi) * exp(-xx^2)
    if (mom==0){
      return(mom0) }
    else if (mom==1){
      return(mean0* mom0 + mom1) }
    else if (mom==2){
      return( (var0+mean0^2)*mom0 + mean0*mom1)}
  }
  
  frmom <- function(mom,meanz,varz,g,PK,momden=FALSE){
    if (identical(momden,FALSE) ){ momden <- mom }
    
    funcresp <- FR$f
    funcrespd <- FR$df
  
    Kint <- PK$table[[as.integer(mom+1)]]
    stdz=sqrt(varz)
    res <- integrate(function(z0){
      (1-funcrespd(z0) *g)^(-momden)*Kint(z0) *dnorm(z0,meanz,stdz)  },-Inf,Inf)$value
    return(res)
  }
  
  meanfield <-function(S,mu,PK){
    v <- PK$avgK/(1+mu*(1-1/S) )
    return(c(1,v,v^2+PK$stdK^2,1 ) )
  }
  
  solve_system <- function(S,mu,sig,gam,PK,correl,maxtrials=100){

    model <- function(x){
      if (identical(FR,FALSE) & identical(PK$dist,FALSE) ){
        gparam <- gaussian_parameters(x,S,mu,sig,gam,PK,correl)
        mean0 <- gparam[1]
        var0 <- gparam[2]
        F1 <- x[1] - erfmom(0,mean0,var0)
        F2 <- x[1] * x[2] - erfmom(1,mean0,var0)
        F3 <- x[1] * x[3] - erfmom(2,mean0,var0)
        F4 <- x[4]*(1-x[4]*gam*sig^2*x[1] ) -1
      }
      else{
        gparam <- gaussian_parameters(x,S,mu,sig,gam,PK,correl)
        mean0 <- gparam[1]
        var0 <- gparam[2]
        zparam <- z_parameters(x,S,mu,sig,gam)
        meanz <- zparam[1]
        varz <- max(0.0001,zparam[2])
        g <- x[1] * sig^2 *gam * x[4]
        g<- min(.9,g)
        F1 <- x[1] - frmom(0,meanz,varz,g,PK)
        F2 <- x[1] * x[2] - frmom(1,meanz,varz,g,PK)
        F3 <- x[1] * x[3] - frmom(2,meanz,varz,g,PK)
        F4 <- x[1] * x[4] - frmom(0,meanz,varz,g,PK,1)
      }
      phimin <-0.001
      cost <- max(0,(phimin-x[1])/phimin ) #max(0,min((0.01-x[1])*10,x[2]^2-x[3]) )
      val <- c(F1 = F1+cost,F2 = F2,F3 = F3, F4=F4) 
      #print(c(x,val, cost ))
      return( val)
    }
    
    ss <- list(root= c( 1,1,1,1),f.root=(100)  )
    failure <- function(s){ max(abs(s$f.root))>0.01  }
    trials=0
    if ((sig>.8 | mu>5.) & trials >0){
      x0 <- solve_system(S,mu/1.1,sig/1.01,gam,PK,correl,maxtrials=2)
    }
    else{
      if (sig>0.1){
        x0 <- solve_system(S,mu/1.1,sig/1.1,0,PK,correl,maxtrials=20)
      }
      else{
        x0 <- meanfield(S,mu,PK)
      }
    }
    while (failure(ss) & trials<maxtrials ){
      try( ss <- multiroot(f = model, start =  x0 * runif(4, .8,1.2), positive=FALSE) ,silent=FALSE)
      trials <- trials +1
      
    }
    if (failure(ss)){
      return( c(0,0,0,0) )
    }
    else{
      return(ss$root)
    }
  }
  
  compute_results <- function(x,S,mu,sig,gam,PK,correl){
    gparam <- gaussian_parameters(x,S,mu,sig,gam,PK,correl)
    u <- gparam[3]
    phi <- x[1]
    N1 <- x[2]
    N2 <- x[3]
    v <- x[4]
    avgN <- N1*phi
    avgN2 <- N2*phi
    stdN <- sqrt(N2-N1^2)
    simpson <- N2/S/phi/N1^2
    zeta <- PK$stdK/PK$avgK
    if(PK$stdK!=0){
      productivity <- (u*avgN2 + mu*avgN^2 +sig^2/PK$stdK^2 
                       *avgN2 * PK$avgK *avgN)/ (1+sig^2/PK$stdK^2*avgN2)
    }
    else{
      productivity <- 1
    }
    press <- 1./( u^2 - phi*sig^2  )
    if(press <0){ press <- Inf }
    variability <- 1/(1- phi* sig^2 * (gam+1)/2 )
    results= c(zeta=zeta,mu=mu,sigma=sig,gamma=gam,
               biomass=avgN*S,sdBiomass=stdN,phi=phi,survivors = S*phi, 
               simpsonD = 1/simpson, productivity = productivity/avgN,
               variability=variability,press=press)
    
    if(correlations){
      results=c(results,cKa=correl[1],sigma_row=correl[2] )
    }
    
    return ( results)
  
  }
  
  ss <- do.call(solve_system,parameters)
  return(do.call(compute_results, c(list(x=ss),parameters)))
  
}

prediction_from_matrix <- function(r,A,D=FALSE,PK=TRUE,correlations=TRUE,FR=FALSE){
  if ( identical(D,FALSE) ){
    D <- diag(A)
    diag(A) <- 0
  }
  return(prediction(measure_parameters(r,A,D,PK,correlations),correlations,FR=FR))
}

from_file <- function(rpath,Apath,Dpath="none"){
  r <- as.matrix(read.table(rpath))
  A <- as.matrix(read.table(Apath))
  if(Dpath=="none"){
    D <- diag(A)}
  else{ 
    D <- as.matrix(read.table(Dpath))
  }
  diag(A) <-0
  return(list(r=r,A=A,D=D))
}


prediction_from_file <- function(rpath,Apath,Dpath="none",PK=TRUE,correlations=TRUE,FR=FALSE)
{
  do.call(prediction_from_matrix,c(from_file(rpath,Apath,Dpath),list(PK=PK,correlations=correlations,FR=FR) ) )
}
