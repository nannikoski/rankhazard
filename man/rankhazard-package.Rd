\name{rankhazard-package}
\alias{rankhazard-package}
\alias{rankhazard}
\docType{package}
\title{Rank-hazard plots}
\description{Rank-hazard plots visualize the relative importance of covariates in a proportional hazards model.  
The key idea is to rank the covariate values and plot the relative hazard as a function of ranks scaled to interval [0,1].  
The relative hazard is the hazard plotted in respect to the reference hazard, which can be e.g. 
the hazard related to the median of the covariate.  
Transformation to scaled ranks allows plotting of covariates measured in different units in the same graph, 
which helps in the interpretation of the epidemiological relevance of the covariates.  
Rank-hazard plots show the difference of hazards between the extremes of the covariate values present 
in the data and can be used as a tool to check if the proportional hazards assumption leads 
to reasonable estimates for individuals with extreme covariate values.  
Alternative covariate definitions or different transformations applied to covariates can be 
also compared using rank-hazard plots.
}
\details{
\tabular{ll}{
Package: \tab rankhazard\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2016-04-01\cr
License: \tab GPL 2\cr
LazyLoad: \tab yes\cr
}
The function \code{rankhazardplot} receives a \code{\link{coxph}} (package \code{survival}) object or a \code{cph} 
(package \code{rms}) object as an argument and creates a rank-hazard plot of the covariates.  
The reference points for the relative hazards and legend texts can be provided as optional arguments.  
Plotting parameters such as, lwd, lty, col and pch are passed to the plotting commands.  
}
\author{
Juha Karvanen <juha.karvanen@iki.fi> , Nanni Koski <nanni.koski@gmail.com> 
}
\references{Karvanen J., Harrell F. E., Jr. 2009 Visualizing covariates in proportional hazards model. 
\emph{Statistics in Medicine}, \bold{28}, 1957--1966.}

\keyword{package}
\keyword{hplot}
\keyword{survival}
\examples{
library(survival)
model <- coxph(Surv(time,status) ~ age + sex + ph.karno, data = cancer, x=TRUE)
rankhazardplot(model, data=cancer)
library(rms)
model2 <- cph(Surv(time,status) ~ age + sex + ph.karno, data = cancer, x=TRUE)
dd <- datadist(cancer)
options(datadist="dd")
rankhazardplot(model2, data=cancer)
}
