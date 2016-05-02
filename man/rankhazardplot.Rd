\name{rankhazardplot}
\alias{rankhazardplot}
\alias{rankhazardplot.cph}
\alias{rankhazardplot.coxph}
\alias{rankhazardplot.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a Rank-hazard plot }
%%%%%%%%%%%%%%%%%%%%%%%%%%%
\description{
Creates a rank-hazard plot.  Plots the relative hazards (or the logarithm of the relative hazards) 
for each covariate of a Cox proportional hazards model fitted by \code{\link{coxph}} or \code{cph}. 
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
rankhazardplot(...)

\method{rankhazardplot}{coxph}(coxphobj, data, select = NULL, refpoints = NULL, 
    CI_level = 0.95, x_CI = NULL, draw.confint = FALSE, legendtext = NULL, 
    axistext = NULL, legendlocation = "top", axistextposition = -0.1, 
    reftick = TRUE, refline = FALSE, refline.col = 1, refline.lwd = 1, 
    refline.lty = 2, ylab = NULL, ylim = NULL, yticks = NULL, 
    yvalues = NULL, plottype = "hazard", na.rm = TRUE, draw = TRUE, 
    return = FALSE, col = NULL, lwd = 1, lty = 1, pch = NULL, axes = TRUE, 
    cex = 1, bg = "transparent", pt.lwd = 1, 
    add = FALSE, graphsbefore = 0, args.legend = NULL, ...)

\method{rankhazardplot}{cph}(cphobj, data, select = NULL, refpoints = NULL, 
    CI_level = 0.95, x_CI = NULL, draw.confint = FALSE, legendtext = NULL, 
    axistext = NULL, legendlocation = "top", axistextposition = -0.1, 
    reftick = TRUE, refline = FALSE, refline.col = 1, refline.lwd = 1, 
    refline.lty = 2, ylab = NULL, ylim = NULL, yticks = NULL, 
    yvalues = NULL, plottype = "hazard", na.rm = TRUE, draw = TRUE, 
    return = FALSE, col = NULL, lwd = 1, lty = 1, pch = NULL, axes = TRUE, 
    cex = 1, bg = "transparent", pt.lwd = 1, 
    add = FALSE, graphsbefore = 0, args.legend = NULL, ...)

\method{rankhazardplot}{default}(x, coefs = NULL, xp = NULL, refvalues = NULL, 
    refpoints = NULL, confinterval = NULL, select = 1, legendtext = NULL, 
    axistext = NULL, legendlocation = "top", axistextposition = -0.1, 
    reftick = TRUE, refline = FALSE, refline.col = 1, refline.lwd = 1, 
    refline.lty = 2, ylab = NULL, ylim = NULL, yticks = NULL, 
    yvalues = NULL, plottype = "hazard", na.rm = TRUE,
    col = NULL, lwd = 1, lty = 1, pch = NULL, axes = TRUE, 
    cex = 1, bg = "transparent", pt.lwd = 1, 
    add = FALSE, graphsbefore = 0, args.legend = NULL, ...)

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{

 \item{coxphobj}{An object of class \code{coxph} created by function \code{\link{coxph}} from the package \code{survival}.  
The object should have been created with \code{x = TRUE}.  Otherwise \code{x_CI} must be given.  }

 \item{cphobj}{An object of class \code{cph} created by function \code{cph} from the package \code{rms}.  
The object should have been created with the option \code{x = TRUE}.  Otherwise \code{x_CI} must be given.  }

 \item{data}{A data frame that contains the covariates in the model.  
It can be the data used in fitting the Cox proportional hazards model or new data.  
It can contain more covariates than there are in the model.  
Used with the argument \code{coxphobj} or \code{cphobj}.  }

 \item{select}{A vector with the order numbers of the covariates to be plotted.  
The order numbers are defined by the order of the covariates in the model.  
It can be used with the argument \code{coxphobj}, \code{cphobj} or \code{confinterval}.  
When used with \code{coxphobj} or \code{cphobj} the default is \code{NULL} and all covariates in the model are plotted. 
With \code{confinterval} the default is the first covariate.  } 

\item{x}{A data frame that contains the covariates to be plotted with one column for each covariate.  
Used with the argument \code{xp} or \code{coefs}.   
When used with \code{xp}, \code{x} defines the values to be printed on the x-axis.  
When used with \code{coefs}, \code{x} the predictions are computed
as a product of \code{coefs} and \code{x}.  
The dimensions of \code{x} and \code{xp} must be the same.  
The number of covariates given in \code{x} and the length of the \code{coefs} must be the same.   
If \code{return = TRUE} \code{x} is returned in the list.  See section \code{Value}.  }

 \item{xp}{A data frame that contains the predictions ("terms") for the covariates to be plotted.  
The order of the covariates must be the same as in \code{x}.  
If \code{return = TRUE} \code{xp} is returned in the list.  See section \code{Value}.  }
 
 \item{coefs}{A vector of Cox regression coefficients for the covariates to be plotted.  
The order of the covariates must be the same as in \code{x}.  }
 
 \item{refpoints}{A vector of reference points given in the same order as the covariates in the model.  
A reference point is the value of the covariate where the reference hazard is calculated to compute 
the relative hazards as a quotient of the hazard and the reference hazard. 
Consequently, at the reference point the value of the relative hazard is 1 (and the value of the logarithm of the relative hazard is 0).   
If the \code{select} argument is in use the reference points are given for the selected covariates only 
and in the same order as the selection made.  
If \code{NULL}, the medians of each covariate are used as reference points.  
With factors the default is the reference level of the factor.  
If the reference point for the selected covariate is \code{NA}, the default is used.  
When plotting the confidence intervals the reference point for factors can be changed only by re-leveling the factor.  }

 \item{refvalues}{A vector of reference values.  
A reference value of a covariate is the predicted value ("terms") at the reference point (see \code{refpoints}).  
Used and needed only with the argument \code{xp}.  
Otherwise calculated by the \code{rankhazardplot} function with the values of the \code{refpoints}.  
If \code{return = TRUE} \code{refvalues} is returned in the list.  See section \code{Value}.  }

 \item{CI_level}{A number between 0 and 1 that defines the level of the confidence interval 
for (the logarithm of) the relative hazard. By default 0.95.} 

 \item{x_CI}{A data frame of the covariate data.  
Needed if the \code{coxphobj$x} or \code{cphobj$x} does not exist.  
The number and the order of the columns must be same as \code{as.data.frame(coxphobj$x)} 
or \code{as.data.frame(cphobj$x)}.  }

\item{col.CI}{Defines the colour of the graph of confidence intervals.  By default the same as colour of the relative hazard.See documentation for \code{\link{par}}.  }

\item{lty.CI}{Defines the type of the graph of confidence intervals.  By default \code{lty + 1}.See documentation for \code{\link{par}}.  }

\item{lwd.CI}{Defines the width of the graph of confidence intervals.  By default the same as width of the relative hazard.See documentation for \code{\link{par}}.  }
 
 \item{draw.confint}{If \code{TRUE} the confidence intervals are plotted.  
Can be used with \code{coxphobj} or \code{cphobj}.  
By default \code{FALSE}.  }

 \item{confinterval}{A list that contains all information to plot confidence intervals.  
If \code{return = TRUE} \code{confinterval} is returned in the list.  See section \code{Value}.  } 
 
 \item{legendtext}{A vector of covariate names for the legend box and the summary.  
If \code{NULL} and \code{axistext = NULL}, the names are from the columns of \code{xp}, \code{coefs}, 
\code{attr(coxphobj$terms,"term.labels")} or \code{attr(cphobj$terms, "term.labels")}.  
If \code{NULL}, \code{axistext} is used if it is given.  
If the \code{select} argument is in use the names are given for the selected covariates only 
and in the same order as the selection made.  }

 \item{axistext}{A vector of covariate names and units for the x-axis.  
If \code{NULL} and \code{legendtext = NULL}, the names are from the columns of \code{x} or the same as in the \code{data}.  
If \code{NULL}, \code{legendtext} is used if it is given.  
If the \code{select} argument is in use the names are given for the selected covariates only 
and in the same order as the selection made.  }

 \item{legendlocation}{A keyword that determines the location where the legend box is printed.  By default "top".
See \code{Details} in the documentation for \code{\link{legend}}.  
If \code{NULL}, the legend box is not printed and it can be added by function \code{\link{legend}}.} 

 \item{axistextposition}{A number that defines the x-coordinate where the axis texts are placed.  
Adjustment is right alignment.  By default -0.1.  
The bigger the size of the plot is, the closer to zero the value can be. 
The maximum is 0.  }

 \item{reftick}{By default \code{TRUE} and the reference tick is emboldened.  
For the relative hazard the tick is at 1 and for the logarithm of the relative hazard at 0.  } 

 \item{refline}{If \code{TRUE} the reference line is drawn.  The line is horizontal 
at the same place as the reference tick (see \code{reftick}).  
By default \code{FALSE}.  } 

 \item{col.refline}{Defines the colour of the reference line, if \code{refline = TRUE}.  By default 1.  
See documentation for \code{\link{par}}.  } 

 \item{lwd.refline}{Defines the width of the reference line, if \code{refline = TRUE}.  By default 1.  
See documentation for \code{\link{par}}.  } 

 \item{lty.refline}{Defines the type of the reference line, if \code{refline = TRUE}.  By default 2.  
See documentation for \code{\link{par}}.  } 

 \item{ylab}{A string that defines the label of the y-axis.  When \code{plottype = "hazard"} the default is "relative hazard".  
When \code{plottype = "loghazard"} the default is the "logarithm of the relative hazard".  See documentation for \code{\link{plot}}.  } 

 \item{ylim}{A numerical vector that defines the range of the y-axis.  
By default the range is the minimum and the maximum of (the logarithm of) the relative hazards.  
See documentation for \code{\link{plot.window}}.  } 

 \item{yticks}{A vector that determines the places for the ticks that are plotted to the y-axis.  
When \code{plottype = "hazard"} it is recommended that evenly calculated ticks are used before and after 
the reference e.g. \code{c(seq(0.5, 1, by=0.1),2:7)}, due to the logarithmic scale of the y-axis.  
If \code{NULL} the ticks are computed using \code{\link{pretty}}.  }  

 \item{yvalues}{A vector that determines which values are printed on the y-axis.  
If \code{NULL} the values of the \code{yticks} are used.  
It is recommended that 1 is in the \code{yvalues} if the \code{plottype = "hazard"} and 0 if the \code{plottype = "loghazard"}.  
In addition to this the \code{yvalues} should be a subset of \code{yticks}.  } 

\item{xtext} {By default \code{TRUE} and the labels and the quartiles of the covariates are printed on the x-axis. If the text on the x-axis should be altered in some way (e.g. \code{font} or \code{cex}) the printing of the labels must be prevented by giving \code{xtext = FALSE}. After that the labels can be added by \code{mtext}}

 \item{axes}{By default \code{TRUE} and the axes are plotted. If \code{FALSE} axes are not plotted. NOTE when adding the axes yourself: On the rank-hazard plot the x-axis is from 0 to 1 and the places for quartiles are 0.25, 0.5 and 0.75. }

 \item{plottype}{A string that defines the scale for the y-axis. 
Either \code{"hazard"} for the relative hazard with log-scale or 
\code{"loghazard"} for the logarithm of the relative hazard with linear scale.} 

 \item{na.rm}{By default \code{TRUE} and only complete cases are plotted.  
Complete cases are cases that have information on all covariates that are used in fitting the model.  
If \code{FALSE} all available cases for each variable are plotted.  } 
 
 \item{draw}{By default \code{TRUE} and a rank-hazard plot and a summary are printed.  
If \code{FALSE}, no output is provided unless \code{return = TRUE}.} 
 
 \item{return}{By default \code{FALSE} and \code{Value} is not returned.  
Used with the argument \code{coxphobj} or \code{cphobj}.  
If \code{TRUE}, \code{x}, \code{xp}, \code{refvalues} and \code{confinterval} are returned as a list.  
See section \code{Value}.  } 
 
 \item{col}{A vector that defines the colours of the lines and the points.  
If the vector is shorter than the number of the covariates to be plotted, the values are repeated.  
See documentation for \code{\link{par}}.  } 

 \item{lwd}{A vector that defines the widths of the lines.  
If the vector is shorter than the number of the covariates to be plotted, the values are repeated.  
See documentation for \code{\link{par}}.  } 

 \item{lty}{A vector that defines the types of the lines.  
If the vector is shorter than the number of the covariates to be plotted, the values are repeated.  
See documentation for \code{\link{par}}.  } 

 \item{pch}{A vector that defines the characters of the points.  
If the vector is shorter than the number of the covariates to be plotted, the values are repeated. If the value is \code{NA} the points are not plotted. See documentation for \code{\link{points}}.} 

 \item{bg}{A vector that defines the fill colour of the point.  Available only for point characters \code{21:25}.
If the vector is shorter than the number of the covariates to be plotted, the values are repeated.  
By default \code{"transparent"}.  
See documentation for \code{\link{points}}.} 

 \item{pt.lwd}{A vector that defines the line width for the drawing symbols.  
If the vector is shorter than the number of the covariates to be plotted, the values are repeated.  
See documentation for \code{lwd} in \code{\link{points}}.} 

 \item{cex}{A vector that defines the size of the points.  
If the vector is shorter than the number of covariates to be plotted, the values are repeated.  
See documentation for \code{\link{plot.default}}.} 

 \item{add}{By default \code{FALSE} and a new plot is created.
If \code{TRUE} graphs are added to the previous plot. In that case the amount of graphs already 
plotted must be specified by the argument \code{graphsbefore}.}
 
 \item{graphsbefore}{A number that specifies how many graphs have already been drawn on the plot. 
If \code{add} = \code{FALSE}, the default is 0.}
  
 \item{args.legend}{A list of additional arguments to be passed to \code{\link{legend}}.
Names of the list are used as argument names. The values are
used when the legend is drawn. For example if \code{legend} and \code{x} are specified in the list, 
those will replace the values of \code{legendtext} and \code{legendlocation} respectively.}

 \item{\dots}{Other arguments to be passed to \code{\link{matplot}} and \code{\link{matlines}} commands.
For example \code{main}. Arguments that are already in use are \code{x}, \code{y}, \code{type}, \code{log}, \code{xlim}, \code{xlab}, \code{xaxt} and \code{yaxt} and using them will give an error.}

}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{

The function \code{rankhazardplot} receives a \code{\link{coxph}} (package \code{survival}) object or a \code{cph} 
(package \code{rms}) object as an argument and creates a rank-hazard plot of the covariates.  
The reference points for the relative hazards and legend texts can be provided as optional arguments.  
Plotting parameters such as, lwd, lty, col and pch are passed to the plotting commands. 

Rank-hazard plots visualize the relative importance of covariates in a proportional hazards model.  
The key idea is to rank the covariate values and plot the relative hazard as a function of ranks scaled to interval [0,1].  
The relative hazard is the hazard plotted in respect to the reference hazard, which can be e.g. 
the hazard related to the median of the covariate. 

The labels on the x-axis show the minimum, the quartiles and the maximum of each covariate.  These are real observations.  
If the quantile would be determined by a mean of two observations, the smaller value is choosed to be the quantile.
However, if the number of the obervations is even, the default reference point is not necessarily a real observation 
as it is calculated as a mean of two middle observations.  Hence, the median shown on the x-axis and the reference point 
can differ even when the default is used. 

Predictions are computed by \code{\link{predict.coxph}}, when the function is called with the argument \code{coxphobj} 
and by \code{predict.rms}, when the function is called with the argument \code{cphobj}.  
Consequently, relative hazards are available for models that are supported by \code{\link{predict.coxph}} or
\code{predict.rms}.  For example the \code{\link{pspline}} transforms are supported by \code{\link{predict.coxph}} 
but not by \code{predict.rms}.

The upper and lower confidence limits of the Cox regression coefficients are used to calculate 
the confidence intervals for the relative hazards.  
Confidence  intervals are only supported for original covariates (same as in the data), 
simply transformed covariates (e.g. \code{log}) and factors. 

Rank-hazard plots can be used to visualize time-dependent models.  In that case plotting can be made using 
\code{coefs} that are Cox regression coefficients of the time-dependent model.  The data matrix \code{x} 
contains the values of covariates at some specific time.  It is also possible to make a rank-hazard plot by the argument
\code{xp}. In that case both \code{x} and \code{xp} must be selected so that they contain information only at some specific time.  
Third way to make the rank-hazard plot is to use the time-dependent model and give the values of the covariates 
at some specific time as the argument \code{data}.  
If the purpose is to compare relative hazards at different times, it is recommended that the same reference points are used in every plot.  
For example the medians of every covariate at the first measurement.


}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{If \code{return = TRUE} the function returns a list that contains all the information needed to draw
 a rank-hazard plot and confidence intervals. The list contains:

\item{x}{A data frame that contains the covariate data.  }

\item{xp}{A data frame that contains the centered predictions for all covariates in the model.  
Calculated by \code{\link{predict.coxph}} or \code{predict.rms}.  }

\item{refvalues}{ A vector that contains the centered predictions that are calculated using the reference point defaults for 
the covariates that don't have a given reference point.  
Calculated by \code{\link{predict.coxph}} or \code{predict.rms}.  }

\item{confinterval}{A list that contains covariate data as a data frame \code{x}, predictions by terms as a data frame \code{xp},
reference values as a vector \code{refvalues}, a lower confidence interval for the predictions as a data frame \code{low}, 
a lower confidence interval for the reference values as a vector \code{lowrefvalues}, an upper confidence interval for the predictions
as a data frame \code{upp} and an upper confidence interval for the reference values as a vector \code{upprefvalues}. 
In the list there is also the vector \code{select_CI} that has indices of the covariates in the model for which the confidence intervals are provided.

The covariates for which the confidence intervals are provided are original (same as in the data), 
simply transformed (e.g. \code{log}) and factors.  The predictions are calculated as 
a product of \code{coefs} and \code{x}.  The upper and lower confidence limits of the Cox regression coefficients 
are used to calculate the confidence intervals for the relative hazards. 
NB: values aren't centered and for that reason e.g. \code{xp} and \code{confinterval$xp} are not the same.  }
}
%%%%%%%%%%%%%%%%%%%%%%%%%

\references{Karvanen J., Harrell F. E., Jr. 2009 Visualizing covariates in proportional hazards model. 
\emph{Statistics in Medicine}, \bold{28}, 1957--1966.}
%%%%%%%%%%%%%%%
\author{Juha Karvanen <juha.karvanen@iki.fi>, Nanni Koski <nanni.koski@gmail.com>}
%%%%%%%%%%%%%%%
\seealso{ \code{\link{coxph}},  \code{cph},\code{\link{predict.coxph}}, \code{predict.rms} }
%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
library(survival)
library(rms)

data(pbc)
# new status variable
pbc$statusbin <- ifelse(pbc$status==0, 0, NA)
pbc$statusbin <- ifelse(pbc$status==2, 1, pbc$statusbin)
pbc$edema2 <- as.factor(pbc$edema)

### different ways to make a rank-hazard plot ###
par(mar = c(4, 5, 4, 2) + 0.1)

coxmodel1 <- coxph(Surv(time, statusbin) ~ age + protime + as.factor(edema),
    data = pbc, x = TRUE)
rankhazardplot(coxmodel1, data = pbc, main = "Rank-hazardplot by coxphobj")

dd <- datadist(pbc)
options(datadist = 'dd')
cmodel1 <- cph(Surv(time, statusbin) ~ age + protime + edema2,
    data = pbc, x = TRUE)
rankhazardplot(cmodel1, data = pbc, main = "Rank-hazardplot by cphobj")

output1 <- rankhazardplot(coxmodel1, data = pbc, draw = FALSE, return = TRUE)
rankhazardplot(x = output1$x, xp = output1$xp, refvalues = output1$refvalues,
    main = "Rank-hazardplot by x and xp")
rankhazardplot(x = output1$x[1:2], coefs = coxmodel1$coef[1:2], 
    main = "Rank-hazardplot by x and coefs \n with returned data")
rankhazardplot(x = pbc[c("age","protime")], coefs = coxmodel1$coef[1:2], 
    main = "Rank-hazardplot by x and coefs \n with original data")

### selecting covariates ###
coxmodel2 <- coxph(Surv(time, statusbin) ~ age + protime + as.factor(edema) +
    bili + albumin + copper + ast + as.factor(stage), data = pbc, x = TRUE)
par(mar = c(9, 5, 4, 2) + 0.1)
rankhazardplot(coxmodel2, data = pbc, main = "Too much information?")
par(mar = c(4, 5, 4, 2) + 0.1)
rankhazardplot(coxmodel2, data = pbc, select = c(1, 4, 5), 
    main = "How to select covariates")

### using graphical parameters ###

# Compare the two following plots
rankhazardplot(coxmodel2, select = c(1, 3, 4, 5), data = pbc, 
    main = "By default")
# X11() # If the plot is too small, legend is not printed in the correct place
par(mar = c(5, 5, 4, 2) + 0.1)
rankhazardplot(coxmodel2, select = c(1, 3, 4, 5), data = pbc,
    legendlocation = "topleft", ylim = c(0.4, 12), 
    ylab = "Relative hazard", yvalues = c(0.4, 1, 2, 4, 6, 10),
    yticks = c(seq(0.4, 1, by = 0.1), 2:10),
    main = "Graphical arguments in use",
    axistext = c("age", "stage", "bilirubin", "albumin"), 
    legendtext = c("age", "factor(stage)", "bilirubin", "albumin"),
    col = c("darkgreen","navyblue","maroon3",1), pch = 18:21, 
    lwd = 2, lty = c(1, 2), cex = 0.9, bg = "yellow", pt.lwd = 2)

# reference line for hazard and loghazard #
par(mfrow = c(1, 2))
rankhazardplot(coxmodel2, select = c(1, 5, 4), data = pbc, refline = TRUE, 
    plottype = "hazard", main = "Reference line at 1")
rankhazardplot(coxmodel2, select = c(1, 5, 4), data = pbc, refline = TRUE,  
    plottype = "loghazard", main = "Reference line at 0")

### comparing covariates from different models ###

# transforms #
# same model as coxmodel2, only bilirubin is transformed by log
logmodel <- coxph(Surv(time, statusbin) ~ age + protime + as.factor(edema) +
    log(bili) + albumin + copper + ast + as.factor(stage), data = pbc, x = TRUE)
# same model as coxmodel2, only a pspline is fitted to bilirubin
coxspline <- coxph(Surv(time, statusbin) ~ age + protime + as.factor(edema) +
    pspline(bili) + albumin + copper + ast + as.factor(stage), data = pbc, x = TRUE)

outputcox <- rankhazardplot(coxmodel2, data = pbc, return = TRUE, draw = FALSE)
outputlog <- rankhazardplot(logmodel, data = pbc, return = TRUE, draw = FALSE)
outputspline <- rankhazardplot(coxspline, data = pbc, return = TRUE, draw = FALSE)

xlog <- data.frame(outputcox$x["bili"], log(outputlog$x["bili"]),
    outputspline$x["bili"])
xp <- data.frame(outputcox$xp["bili"], outputlog$xp["log(bili)"],
    outputspline$xp["pspline(bili)"])
ref <- c(outputcox$ref["bili"], outputlog$ref["log(bili)"],
    outputspline$ref["pspline(bili)"])
par(mfrow = c(1, 1))
par(mar = c(3, 5, 4, 2) + 0.1)
rankhazardplot(x = xlog, xp = xp, refvalues = ref, 
    legendtext = c( "bili","log(bili)", "pspline(bili)"), 
    axistext = c("bili", "log(bili)", "bili"),
    main = "Transforming has a great impact on \n interpreting the effect of the bilirubin")

# factors #
# same model as coxmodel2, only age is left out
# how does it affect relative hazards for stage?
coxmodel3 <- coxph(Surv(time, statusbin) ~ protime + as.factor(edema) +
    log(bili) + albumin + copper + ast + as.factor(stage), data = pbc, x = TRUE)
outputcox2 <- rankhazardplot(coxmodel3, data = pbc, return = TRUE, draw = FALSE)
xp <- data.frame(outputcox$xp["as.factor(stage)"],outputcox2$xp["as.factor(stage)"])
x <- data.frame(outputcox$x["stage"],outputcox2$x["stage"])
ref <- c(outputcox$ref["as.factor(stage)"],outputcox2$ref["as.factor(stage)"])

par(mar = c(2, 5, 4, 2) + 0.1)
rankhazardplot(x = x, xp = xp, refvalues = ref, legendlocation = "bottomright",
    legendtext = c("age in the model","age not in the model"), 
    axistext = c("stage","stage"),
    main = "Comparing the relative risk of 'stage' between
    similar models that differ by covariate 'age'")

### changing reference points ###

# with no model object #
coefs <- c(coxmodel2$coef["bili"],logmodel$coef["log(bili)"])
par(mar = c(3, 5, 4, 2) + 0.1)
rankhazardplot(x = xlog[1:2], coefs = coefs[1:2], refpoints = c(1.2, log(1.2)), 
    legendtext = c("bili", "log(bili)"), 
    main = "Reference point is 1.2, the maximum \n value for 'normal' total bilirubin")

# factors with non-numerical levels#
coxmodel4 <- coxph(Surv(time, statusbin)~age + sex, data = pbc, x = TRUE)
par(mfrow = c(1, 2))
rankhazardplot(coxmodel4, data = pbc, main = "Reference points by default", 
    refline = TRUE, ylim = c(0.4, 3.8))
rankhazardplot(coxmodel4, refpoints = c(40, "f"), data = pbc,
    main = "Different reference points", refline = TRUE,ylim = c(0.4, 3.8))

par(mfrow = c(1, 1))
# with select argument and changing only part of reference points #
rankhazardplot(coxmodel2, data = pbc, select = c(7, 1), refpoints = c(100, NA),
    ylim = c(0.5, 4.3), refline = TRUE, main = "Reference point for age by default")

# using xp argument #
output2 <- rankhazardplot(coxmodel1, refpoints = c(40, 10,0) , data = pbc, 
    draw = FALSE, return = TRUE)
rankhazardplot(x = output2$x, xp = output2$xp, refvalues = output2$refvalues, 
    main = "How to change the reference \n points when using xp")

### confidence intervals ###
par(mar = c(2, 5, 4, 2) + 0.1)
rankhazardplot(confinterval = output1$conf, 
    main = "By argument confinterval, \n 95 per cent confidence intervals")
rankhazardplot(coxmodel2, data = pbc, confint = TRUE, select = 1,
    refpoint = 40, main = "By argument confint and \n changing reference point")

### data in start-stop format ###
data(cgd)
timemodel <- coxph(Surv(tstart, tstop, status) ~ treat + height + steroids, data = cgd, x = TRUE)
# steroids and height are in the model only to make the example plot more interesting
rankhazardplot(timemodel, data = cgd[cgd$enum == 1,], main = "Covariate values at study entry")


%function(...) UseMethod("rankhazardplot")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{survival}
