rankhazardplot.coxph <- function (
    coxphobj, data = NULL, select = NULL, refpoints = NULL, 
    CI_level = 0.95, x_CI = NULL, draw.confint = FALSE, legendtext = NULL, 
    axistext = NULL, legendlocation = "top", axistextposition = -0.1, 
    reftick = TRUE, refline = FALSE, refline.col = 1, refline.lwd = 1, 
    refline.lty = 2, ylab = NULL, ylim = NULL, yticks = NULL, 
    yvalues = NULL, plottype = "hazard", na.rm = TRUE, draw = TRUE, 
    return = FALSE, col = NULL, lwd = 1, lty = 1, pch = NULL, axes = TRUE,
    cex = 1, bg = "transparent", pt.lwd = 1, add = FALSE, graphsbefore = 0, args.legend = NULL, ...)		
{
  if (!is.data.frame(data)) 
    stop("Covariate data need to be provided as a data frame by the argument 'data'.")

  if (!is.matrix(coxphobj$x) & is.null(x_CI))
    stop("To calculate confidence intevals covariate data need to be provided either as argument 'x_CI' or as 'coxphobj$x'.")
  if (is.null(x_CI)) x_CI <- as.data.frame(coxphobj$x)
  if (!is.data.frame(x_CI)) stop("'x_CI' must be a data frame.")
  
  term_labels <- attr(coxphobj$terms, "term.labels")

  if (is.null(select))
    select <- 1:length(term_labels)
  if (max(select) > length(term_labels))
    stop("There are fewer covariates in the model than selected.")
    

  trans_var <- which(!is.element(term_labels, names(data)))	# checks which term labels are different as in the data
  # changes the labels of the transformed variables to get the same as in the data
  labels <- sapply(strsplit(term_labels[trans_var], ',', fixed = TRUE), '[[', 1)	# extracts first string before ','
  labels <- gsub("^.*\\(", "", labels) 		# deletes all characters before '('
  labels <- gsub(")", "", labels, fixed = TRUE) 				# deletes all ')'
  # combines original and changed labels
  data_labels <- term_labels						
  data_labels[trans_var] <- labels
  x <- data[data_labels] # yritä rikkoa#########################

  if (na.rm) x <- na.omit(x)

  factorlevels <- coxphobj$xlevels
  factorlabs <- names(factorlevels)
  factors <- which(is.element(term_labels, factorlabs))
  nonfactors <- which(!is.element(term_labels, factorlabs))
#browser()
  refs <- data.frame(setNames(replicate(length(data_labels), numeric(0), simplify = F), data_labels))

  for(i in nonfactors)
    refs[1, i] <- median(x[, i], na.rm = TRUE)

  j <- 1
  for(i in factors){
    refs[1, i] <- factorlevels[[j]][1]
    refs[i] <- factor(refs[i], levels = factorlevels[[j]])
    j <- j + 1
  }
  
     if(!is.null(refpoints)){
     if(length(refpoints) != length(select)) 
       stop("The length of 'refpoints' must be the same as the number of covariates to be plotted.")
     change <- which(!is.na(refpoints))
     # testi että ei-faktoreilla on numeeriset refpointit
#browser()
     if (any(is.na(as.numeric(refpoints[intersect(select[change], nonfactors)])))) warning("'refpoints' must be numeric for variables which are not factors")
# testi että faktorien refpointit kuuluu faktorin leveleihin
     j <- 1
     for (i in factors){
       if (is.element(i, select[change]))
         if (!is.element(refpoints[which(select == i)], factorlevels[[j]])) 
           warning(paste("The given value in 'refpoints' for factor", factorlabs[j], "is not one of the levels of the factor."))
       j <- j + 1
     }  
     if(is.numeric(refpoints)){ #testi että ei ole faktoreita?
        refs[1, select[change]] <- refpoints[change]
     }else{
        refs[1, intersect(factors, select[change])] <- refpoints[is.element(select, intersect(factors, select[change]))] #kyllä tämä toimii
        refs[1, intersect(nonfactors, select[change])] <- as.numeric(refpoints[is.element(select, intersect(nonfactors, select[change]))])
      }
}

   predictions <- predict(coxphobj, type = "terms", newdata = x)
   pred_refvalues <- predict(coxphobj, type = "terms", newdata = refs)

   refvalues <- as.vector(pred_refvalues)
   names(refvalues) <- attr(pred_refvalues, "dimnames")[[2]]

   xp <- as.data.frame(predictions)

### Calculating the confidence intervals ###
    
   refs[factors] <- 0
   refs <- as.vector(as.matrix(refs))

   if (CI_level <= 0 || 1 <= CI_level) stop ("'CI_level' must be a number between 0 and 1.")
   confinterval <- rankhazard_CI.coxph(coxphobj, x_CI, refs, CI_level)

   select_CI <- confinterval$select_CI
   selecttext <- select

   if (draw.confint){
     CI <- confinterval
     select <- which(is.element(select_CI, select))
     selecttext <- select_CI[select]
     if(length(select) == 0)
       stop("Confidence intevals cannot be caluclated for selected covariates.")
    } else {
        CI<- NULL
    }

    if (!is.null(legendtext) & is.null(axistext))
      axistext <- legendtext
    if (is.null(legendtext) & !is.null(axistext))
      legendtext <- axistext
    
    if (is.null(legendtext)){
      legendtext <- term_labels[selecttext]
      axistext <- data_labels[selecttext]			
    }

    if (draw)
      rankhazardplot.default( 						
        x = x[select], xp = xp[select], refvalues = refvalues[select], 					
        legendtext = legendtext, axistext = axistext,
        na.rm = na.rm, select = select, confinterval = CI, 
        legendlocation = legendlocation, axistextposition = axistextposition, 
        reftick = reftick, refline = refline, refline.col = refline.col, 
        refline.lwd = refline.lwd, refline.lty = refline.lty, ylab = ylab, 
        ylim = ylim, yticks = yticks, yvalues = yvalues, plottype = plottype, 
        col = col, lwd = lwd, lty = lty, pch = pch, axes = axes,
        cex = cex, bg = bg, pt.lwd = pt.lwd, add = add, graphsbefore = graphsbefore, 
        args.legend = args.legend, ...)

    if(return)
      return(list(x = x, xp = xp, refvalues = refvalues, confinterval = confinterval))
}