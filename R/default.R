rankhazardplot <- function(...) UseMethod("rankhazardplot")


rankhazardplot.default <- function (
    x, coefs = NULL, xp = NULL, refvalues = NULL, refpoints = NULL,
    confinterval = NULL, select = 1, legendtext = NULL, 
    axistext = NULL, legendlocation = "top", axistextposition = -0.1, 
    reftick = TRUE, refline = FALSE, refline.col = 1, refline.lwd = 1, 
    refline.lty = 2, ylab = NULL, ylim = NULL, yticks = NULL, 
    yvalues = NULL, plottype = "hazard", na.rm = TRUE,
    col = NULL, lwd = 1, lty = 1, pch = NULL, 
    cex = 1, bg = "transparent", pt.lwd = 1, add = FALSE, graphsbefore = 0, ...)				
{
    if(!is.null(confinterval)){
        x <- confinterval$x
        if (na.rm) x <- na.omit(x)
        x <- confinterval$x[select]
        xp <- confinterval$xp
        if (na.rm) xp <- na.omit(xp)
        xp <- confinterval$xp[select]
        refvalues <- confinterval$refvalues[select]
    }

    if (na.rm) x <- na.omit(x)
    if (na.rm & !is.null(xp)) xp <- na.omit(xp)

    n <- dim(x)[1]	# number of observations
    m <- dim(x)[2]	# number of covariates

    if (!identical(plottype, "hazard") & !identical(plottype, "loghazard")) 		
        stop("Unknown plottype")
    if (is.null(xp) & is.null(coefs)) 					
        stop("Either coefs or xp must be provided.")		
    if (is.null(refvalues) & !is.null(xp)) 				
        stop("When xp is given, also refvalues are required.")	

    if(is.null(refvalues) & is.null(refpoints)){	
        refpoints <- apply(x, 2, median, na.rm = TRUE)	
        refvalues <- coefs*refpoints
    }

    if(is.null(refvalues) & !is.null(refpoints))	
        refvalues <- coefs*refpoints

    if (is.null(xp)) 
        xp <- as.data.frame(t(coefs * t(x))) 

    lwd <- rep(lwd, length.out = m)
    lty <- rep(lty, length.out = m)
    cex <- rep(cex, length.out = m)
    bg <- rep(bg, length.out = m)
    pt.lwd <- rep(pt.lwd, length.out = m)
    
    if (!add) graphsbefore = 0 #makes sure that 'graphsbefore' is only in use with 'add = TRUE'

    if (is.null(pch)){pch <- seq(0, m - 1) + graphsbefore} 		
    else{pch <- rep(pch, length.out = m)}							
    if (is.null(col)) {col <- 1:m + graphsbefore }				
    else{ col <- rep(col, length.out = m)	}	
    
    if (is.null(legendtext) & !is.null(axistext)) 
        legendtext <- axistext	
    if (!is.null(legendtext) & is.null(axistext)) 
        axistext	<- legendtext	 
    if (is.null(legendtext) & is.null(axistext) & !is.null(names(xp))) 
        legendtext <- names(xp)			
    if (is.null(legendtext) & is.null(axistext) & !is.null(names(coefs))) 
        legendtext <- names(coefs)
    if (is.null(axistext) & !is.null(colnames(x))) 
        axistext <- colnames(x)			
    if (is.null(axistext)) 
        axistext <- legendtext

    ones <- matrix(1, nrow = n, ncol = 1)	
    y <- xp - ones %*% refvalues

    if(!is.null(confinterval)){
        upp_ci <- confinterval$upp - ones %*% confinterval$upprefvalues
        upp_ci <- upp_ci[select]
        low_ci <- confinterval$low - ones %*% confinterval$lowrefvalues
        low_ci <- low_ci[select]
    }

    if (identical(plottype, "hazard")){
        y <- exp(y)
    }
    if (identical(plottype, "hazard") & !is.null(confinterval)){
        low_ci <- exp(low_ci) 
        upp_ci <- exp(upp_ci)
    }
    yrange <- y    # makes sure that confidence intervals fit to the screen
    if(!is.null(confinterval)){
        yrange <- as.data.frame(c(y, low_ci, upp_ci))	 
    }  

    if (length(ylim)!= 2){
        maxy <- max(yrange, na.rm = TRUE)
        miny <- min(yrange, na.rm = TRUE) 
    }else{
        maxy <- ylim[2]
        miny <- ylim[1]
    }

    if (identical(plottype, "hazard")) {	
        if(is.null(ylab)) ylab <- "relative hazard"					
        if (is.null(yticks))
            yticks <- c(pretty(c(miny, 1)), pretty(c(1, maxy)))	
        reftickvalue <- 1
        logvar = "y"
    }    											
    if (identical(plottype, "loghazard")) {
        if(is.null(ylab)) ylab <- "logarithm of the relative hazard"										
        if (is.null(yticks))
            yticks <- c(pretty(c(miny, 0)), pretty(c(0, maxy)))
        reftickvalue <- 0	
        logvar = ""
    }

    if (is.null(yvalues)) yvalues <- yticks
    quantiles <- c(0, 0.25, 0.5, 0.75, 1)	
    orders <- apply(x, 2, order) # orders = ind (later)

###lisätty alkaa###
    scaleranks <- x
    y_ord <- y
    rank_quantile <- matrix(ncol=m, nrow=5)
    y_points <- matrix(ncol=m, nrow=5)
    
    na_sum <- colSums(is.na(x))
 
    for(i in 1:m){
      scaleranks[i] <- c(seq(0, 1, length = n - na_sum[i]), rep(NA, na_sum[i]))
      y_ord[i] <- y[orders[,i],i]
      rank_quantile[,i] <-quantile(1:(n - na_sum[i]), probs = quantiles)
      y_points[,i] <- y_ord[,i][rank_quantile[,i]]
    }

    matplot(scaleranks, y_ord, type="l", log=logvar, ylim=c(miny, maxy), xlim=c(0,1), ylab=ylab, xlab="", xaxt="n", yaxt = "n",
            col=col, lty=lty, lwd=lwd, add = add, ...)
    matpoints(quantiles, y_points, pch=pch, col=col, cex=cex, bg=bg, lwd=pt.lwd)
    
    if (!is.null(confinterval)){
      low_ci_ord <- low_ci
      upp_ci_ord <- upp_ci
      low_ci_points <- matrix(ncol=m, nrow=5)
      upp_ci_points <- matrix(ncol=m, nrow=5)
      
      for(i in 1:m){
        low_ci_ord[i] <- low_ci[orders[ , i], i]
        upp_ci_ord[i] <- upp_ci[orders[ , i], i]
        low_ci_points[ , i] <- low_ci_ord[ , i][rank_quantile[ , i]]
        upp_ci_points[ , i] <- upp_ci_ord[ , i][rank_quantile[ , i]]
      }
     
      matlines(scaleranks, low_ci_ord, type = "l",col = col, lty = lty + 1, lwd = lwd, add = add, ...) 
      matlines(scaleranks, upp_ci_ord, type = "l",col = col, lty = lty + 1, lwd = lwd, add = add, ...) 
      matpoints(quantiles, low_ci_points, pch = pch, col = col, cex = cex, bg = bg, lwd = pt.lwd)
      matpoints(quantiles, upp_ci_points, pch = pch, col = col, cex = cex, bg = bg, lwd = pt.lwd)
      
    }
    
    for(i in 1:m){
      xlabels <- x[orders[rank_quantile[,i],i], i]    # quantiles for covariate i
      if (is.numeric(xlabels)) xlabels <- signif(xlabels, 3) #rounds numeric labels
      mtext(side = 1, at = c(axistextposition, quantiles),    
            adj = c(1,rep(0.5, length(quantiles))), text = c(axistext[i], as.character(xlabels)), line = i + graphsbefore)
    }
    
    if (!add){
      axis(1, at = quantiles, labels = FALSE)    # marks ticks on x-axis
      axis(2, at = yticks, labels = FALSE)    # marks ticks on y-axis
      axis(2, at = yvalues, labels = as.character(yvalues))    # marks values on y-axis
      
      if (reftick)    # eboldens the reference tick
        axis(2, at = reftickvalue, labels = FALSE, lwd.ticks = 2)
      
      if (refline)    # draws the reference line
        abline(h = reftickvalue,  col = refline.col, lty = refline.lty, lwd = refline.lwd)
      
      legend(legendlocation, legend = legendtext, col = col, lwd = lwd, 
             pch = pch, lty = lty, bty = "n", pt.cex = cex, pt.lwd = pt.lwd, pt.bg = bg)
    }

    ### Output to console ####
    A <-matrix(0, m, 5)
    colnames(A) <- c("Min.", "1st Qu.", "Median" , "3rd Qu.", "Max.")
    rownames(A) <- legendtext
    
    for(i in 1:m)
      A[i,] <- quantile(y[, i], probs = quantiles, na.rm = TRUE) #osaako median ottaa quantiilit
    
    cat("Y-axis range: ", signif(c(miny, maxy), 3), "\n", "\n")
    if (identical(plottype, "hazard")) cat("Relative hazards for each covarite:", "\n")
    if (identical(plottype, "loghazard")) cat("Logarithm of the relative hazards for each covarite:", "\n")
    print(signif(A, 3))

    if (!is.null(confinterval)){

        cat("\n")
        if (identical(plottype, "hazard")) 
            cat("Relative hazards for the confidence intervals of each covariate:", "\n")
        if (identical(plottype, "loghazard")) 
            cat("Logarithm of the relative hazards for the confidence intervals of each covariate:", "\n")
        B <- matrix(0, 2 * m, 5)
        colnames(B) <- c("Min.", "1st Qu.", "Median" , "3rd Qu.", "Max.")
	  low_legend <- paste("Low", legendtext, sep = "_")
	  upp_legend <- paste("Upp", legendtext, sep = "_")
        rownames(B)[2 * 1:m] <- upp_legend
        rownames(B)[2 * 1:m - 1] <- low_legend

        for(i in 1:m){
            B[2 * i - 1,] <- quantile(low_ci[, i], probs = quantiles, na.rm = TRUE)
            B[2 * i,] <- quantile(upp_ci[, i], probs = quantiles, na.rm = TRUE)
        }
    
        print(signif(B, 3))
###############
    }		


}

