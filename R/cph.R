#5.2.2016

rankhazardplot.cph <- function (
    cphobj, data = NULL, select = NULL, refpoints = NULL, 
    CI_level = 0.95, x_CI = NULL, confint = FALSE, legendtext = NULL, 
    axistext = NULL, legendlocation = "top", axistextposition = -0.1, 
    reftick = TRUE, refline = FALSE, refline.col = 1, refline.lwd = 1, 
    refline.lty = 2, ylab = NULL, ylim = NULL, yticks = NULL, 
    yvalues = NULL, plottype = "hazard", na.rm = TRUE, draw = TRUE, 
    return = FALSE, col = NULL, lwd = 1, lty = 1, pch = NULL, 
    cex = 1, bg = "transparent", pt.lwd = 1, ...)					
{
    if (is.null(data)) 
        stop("Covariate data need to be provided as argument data.")

    if (is.null(cphobj$x) & is.null(x_CI))
        stop("To calculate confidence intevals covariate data need to be provided either as argument x_CI or as cphobj$x.")

    if (is.null(x_CI))
        x_CI <- as.data.frame(cphobj$x)
 
    term_labels <- cphobj$Design$name 

    if (is.null(select))
        select <- 1:length(term_labels)

    data_labels <- term_labels						

    x <- data[data_labels]

    if (na.rm) x <- na.omit(x)

    factorlevels <- cphobj$Design$parms
    factorlabs <- names(factorlevels)
    factors <- which(is.element(term_labels, factorlabs))
    nonfactors <- which(!is.element(term_labels, factorlabs))

    refs <- data.frame(setNames(replicate(length(data_labels), numeric(0), simplify = F), data_labels))

    for(i in nonfactors)
    	refs[1, i] <- median(x[, i], na.rm = TRUE)

    j <- 1
    for(i in factors){
        refs[1, i] <- factorlevels[[j]][1]
        refs[i] <- factor(refs[i], levels  = factorlevels[[j]])
        j <- j + 1
    }

    if(!is.null(refpoints)){
        change <- which(!is.na(refpoints))
        if(is.numeric(refpoints)){
            refs[1, select[change]] <- refpoints[change]
        }else{
            refs[1, intersect(factors, select[change])] <- refpoints[is.element(select, intersect(factors, select[change]))]
            refs[1, intersect(nonfactors, select[change])] <- as.numeric(refpoints[is.element(select, intersect(nonfactors, select[change]))])
        }
    }

    predictions <- predict(cphobj, type = "terms", newdata = x)
    pred_refvalues <- predict(cphobj, type = "terms", newdata = refs)

    refvalues <- as.vector(pred_refvalues)
    names(refvalues) <- attr(pred_refvalues, "dimnames")[[2]]

    xp <- as.data.frame(predictions)

### Calculating the confidence intervals ###

    coefslow <- confint(cphobj, level = CI_level)[, 1]
    coefsupp <- confint(cphobj, level = CI_level)[, 2]

    refs[factors] <- 0
    refs <- as.vector(as.matrix(refs))

    Values <- cph_CI(cphobj, x_CI, cphobj$coef, refs)
    CIlow <- cph_CI(cphobj, x_CI, coefslow, refs)
    CIupp <- cph_CI(cphobj, x_CI, coefsupp, refs)

    confinterval <- list(x = Values$x, xp = Values$xp, refvalues = Values$refvalues, low = CIlow$xp, lowrefvalues = CIlow$refvalues, upp = CIupp$xp, upprefvalues = CIupp$refvalues)
    select_CI <- Values$select_CI
    selecttext <- select

    if (confint){
        CI <- confinterval
        select <- which(is.element(select_CI, select))
        selecttext <- select_CI[select]
        if(length(select) == 0)
            stop("Confidence intevals cannot be caluclated for selected covariates.")
    } else {
        CI <- NULL
    }

    if (!is.null(legendtext) & is.null(axistext)){
    axistext <- legendtext
    }
    if (is.null(legendtext) & !is.null(axistext)){
    legendtext <- axistext
    } 
    if (is.null(legendtext)){
        legendtext <- attr(cphobj$terms, "term.labels")[selecttext]
        axistext <- term_labels[selecttext]		
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
            col = col, lwd = lwd, lty = lty, pch = pch, 
            cex = cex, bg = bg, pt.lwd = pt.lwd, ...)

    if (return)
        return(list(x = x, xp = xp, refvalues = refvalues, confinterval = confinterval))
}
