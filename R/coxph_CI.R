#5.2.2016

coxph_CI <- function(coxphobj, x, coefs, refpoints){

# This function calculates predictions as a product of coefs and x, and
# reference values as a product of coefs and refvalues.
# Factors are given in a dummy format but returned as one variable.
# The function returns data as x, predictions as xp, reference values as refvalues
# and indices of covariates for which the predictions can be calculated by this function as
# select_CI.

    factorlevels <- coxphobj$xlevels
    factorlabs <- names(factorlevels)
    xp <- as.data.frame(t(coefs * t(x)))

    covariatelabs <- attr(coxphobj$terms, "term.labels")	
    factors <- which(is.element(covariatelabs, factorlabs))	
    columns <- sapply(coxphobj$assign, length)            
    orig_var <- which(columns == 1)
    select <- sort(union(orig_var, factors))        #for these covariates confidence intervals can be calculated
    indices <- sapply(coxphobj$assign, "[[", 1)        #indices where values for each covariate are/begin

    xp[factorlabs] <- 0        #these covariates don't exist yet because factors are in a dummy format
    x[factorlabs] <- 0

    j <- 1
    for(i in factors){    # levels of factors are combined to one variable

        if (columns[i] > 1)
            xp[names(columns)[i]] <- apply(xp[, indices[i] + 0:(columns[i] - 1)], 1, sum) 
        else xp[names(columns)[i]] <- xp[, indices[i]]

        for (l in 1:(columns[i])){   #the values of different levels are copied into same variable
            x[factorlabs[j]] <- ifelse(x[, indices[i] + l - 1] == 1, factorlevels[[j]][l + 1], x[, factorlabs[j]])
        }
        #the cases with xp == 0 have the value of the reference level
        x[factorlabs[j]] <- ifelse(xp[, factorlabs[j]] == 0, factorlevels[[j]][1], x[, factorlabs[j]]) 
        j <- j + 1
    } 

    #the factors in the model are coerced as factors in x

    j <- 1
    for(i in factors) {
        xfactor <- as.factor(x[, factorlabs[j]])
        x[factorlabs[j]] <- relevel(xfactor, ref = factorlevels[[j]][1])
        j <- j + 1
    }

    covariatelabs <- covariatelabs[select]

    xp <- xp[covariatelabs] 		
    x <- x[covariatelabs]			

    refvalues <- coefs[indices[select]] * refpoints[select]

    names(refvalues) <- covariatelabs			
    refvalues[factorlabs] <- 0          # the reference value for factors is zero as coefficient for the reference level is zero

    return(list(x = x, xp = xp, refvalues = refvalues, select_CI = select))
}

