library(rootSolve)
source("./comparison.R")

# This code is provided as part of the Supporting Information for the article
#     "Generic assembly patterns in large ecological communities" (2017)
#     M. Barbier, J.-F. Arnoldi, G. Bunin and M. Loreau.


# Basic example of prediction from example files
prediction_from_file("r.txt","A.txt","D.txt",PK=FALSE)

# Prediction example with only the four basic reference parameters (not the additional correlations)
prediction_from_file("r.txt","A.txt","D.txt",correlations=FALSE,PK=FALSE )

# Prediction example with explicit functional response
# FR should be a list containing the functional response f(x) and its derivative f'(x)
FR <- list(f=function(x){x},df=function(x){1})   # Linear functional response
prediction_from_file("r.txt","A.txt","D.txt",FR=FR )

# Basic example of comparison between data and predictions
comparison_from_file("B.txt","r.txt","A.txt","D.txt",PK=TRUE)

