library(rootSolve)
source("./prediction.R")

# This code is provided as part of the Supporting Information for the article
#     "Generic assembly patterns in large ecological communities" (2017)
#     M. Barbier, J.-F. Arnoldi, G. Bunin and M. Loreau.

# The function "comparison_from_file(Bpath,rpath,Apath,Dpath) takes paths to
# the final abundances B, and the coefficients r, A and D.
# It returns a list
#     Comparison: dataframe containing predictions and simulation outcomes
#         for various community properties (see "prediction" for a list)
#     Error: average error between predictions and simulation outcomes

# NB: For now, variability is not included in the comparison
# because we could not find an R implementation of the algorithm
# for solving the continuous-time Lyapunov equation.

compare <- function(predictions,results){
  nm <- names(results)
  errors <-abs(predictions[nm]/results[nm]-1)
  dataf <- data.frame(predictions[nm],results,error=errors)
  
  return(dataf)
}

measure_results <- function(B,r,A,D,death=10^-8){
  K <-r/D
  alpha <- A / rep(D, each = nrow(A))
  
  S <- length(B)
  alive <- (B>death)
  survivors <- sum(alive)
  phi <- survivors *1. / S
  N1 <- mean(B[alive])
  N2 <- mean( (B[alive])^2  )
  stdN <- sqrt(N2-N1^2)
  simpson <- N2/S/phi/N1^2
  productivity <- mean(B*K)/mean(B)
  
  J <-  (-alpha-diag( rep(1,S)  ))[alive,alive]
  invJ <-  solve(J) 
  press <-norm(invJ,type="F")^2/survivors
  
  results <- c(biomass=N1*S*phi, sdBiomass=stdN, phi=phi,survivors=survivors,
               simpsonD=1/simpson,productivity=productivity, press=press)
  return(results )
  
}

comparison_from_matrix <- function(B,r,A,D,PK=TRUE,correlations=TRUE,FR=FALSE){
  predictions <- prediction_from_matrix(r,A,D,PK,correlations,FR)
  measures <- measure_results(B,r,A,D)
  compare(predictions,measures)
  
}

comparison_from_file <- function(Bpath,rpath,Apath,Dpath,PK=TRUE,correlations=TRUE,FR=FALSE){
  mat <- from_file(rpath,Apath,Dpath)
  B <- as.matrix(read.table(Bpath))
  comparison_from_matrix(B,mat$r,mat$A,mat$D,PK,correlations,FR)

}

