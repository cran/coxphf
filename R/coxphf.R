".First.lib" <-
function(lib, pkg) 
		library.dynam("coxphf", pkg, lib)

"coxphf" <-
function(
 formula=attr(data, "formula"),	# formula where .time. represents time for time-interactions
 data=sys.parent(),
 pl=TRUE,
 alpha=0.05,
 maxit=50,
 maxhs=5,
 epsilon=1e-6,
 maxstep=2.5,
 firth=TRUE
){
### by MP und GH, 2006
### Surv-Objekt entweder usual  (z.B. Surv(time,event)~A+B+C+D) oder
### Counting-Process-Style!! (z.B. Surv(start,stop,event)~A+B+C+D)
### event: 1=dead/event, 0=censored
###
        n <- nrow(data)

	  obj <- decomposeSurv(formula, data, sort=TRUE)
        NTDE <- obj$NTDE
	  mmm <- cbind(obj$mm1, obj$timedata)
        
	  cov.name <- obj$covnames

        k <- ncol(obj$mm1)          # number covariates
        ones <- matrix(1, n, k+NTDE)
          
	  ## standardisierung
	  sd1 <- sd(obj$mm1)
	  sd2 <- sd(obj$timedata)
	  Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
	  obj$mm1 <- scale(obj$mm1, FALSE, sd1)
	  obj$timedata <- scale(obj$timedata, FALSE, sd2)
	  mmm <- cbind(obj$mm1, obj$timedata)
	   
	  CARDS <- cbind(obj$mm1, obj$resp, ones, obj$timedata)	  
        PARMS <- c(n, k, firth, maxit, maxhs, maxstep, epsilon, 0, 0, 0, 0, 0, 0, NTDE)
        IOARRAY <- rbind(rep(1, k+NTDE), matrix(0, 2+k+NTDE, k + NTDE))
        if(NTDE>0)
                IOARRAY[4, (k+1):(k+NTDE)] <- obj$timeind
        storage.mode(CARDS) <- "double"
        storage.mode(PARMS) <- "double"
        storage.mode(IOARRAY) <- "double"

        ## --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
        value <- .Fortran("firthcox",
                CARDS,
                outpar = PARMS,
                outtab = IOARRAY, PACKAGE=coxphf)
        if(value$outpar[8])
                warning("Error in routine FIRTHCOX; parms8 <> 0")
        outtab <- matrix(value$outtab, nrow=3+k+NTDE) #

        # --------------- Verarbeitung des Outputs zu Objekt <fit> der Klasse coxphf -------
        coef.orig <- outtab[3,  ]
        coefs <- coef.orig / Z.sd
        covs <- matrix(outtab[4:(k+3+NTDE), ], ncol=k+NTDE) / (Z.sd %*% t(Z.sd))

        dimnames(covs) <- list(cov.name, cov.name)
        vars <- diag(covs)
        names(coefs) <- cov.name
        fit <- list(coefficients = coefs, alpha = alpha, var = covs, df = k+NTDE,
                        loglik = value$outpar[12:11], iter = value$outpar[10],
                        method.ties = "breslow", n = n, ##terms = terms(formula),
                        y = obj$resp, formula = formula, call = match.call())
        fit$means <- apply(mmm, 2, mean)
        fit$linear.predictors <- as.vector(scale(mmm, fit$means, scale=FALSE) %*% coefs)
        if(firth)
                fit$method <- "Penalized ML"
          else
                fit$method <- "Standard ML"


        # --------------- Aufruf Fortran - Makro PLCOMP ------------------------------------
        if(pl) {
                PARMS <- c(PARMS[1:7], qchisq(1-alpha, 1), 0, 0, 0, 0, 0, NTDE)
                IOARRAY <- rbind(rep(1, k+NTDE), rep(0, k+NTDE), coef.orig, matrix(0, 5, k+NTDE))
                if(NTDE>0)
                		IOARRAY[4,(k+1):(k+NTDE)] <- obj$timeind
                storage.mode(PARMS) <- "double"
                storage.mode(IOARRAY) <- "double"

                value <- .Fortran("plcomp",
                        CARDS,
                        outpar = PARMS,
                        outtab = IOARRAY, PACKAGE=coxphf)
                if(value$outpar[9])
                        warning("Error in routine PLCOMP: parms9 <> 0")
                outtab <- matrix(value$outtab, nrow = 8, ncol = k+NTDE)
                fit$method.ci <- "Profile Likelihood"
                fit$ci.lower <- exp(outtab[4,  ] / Z.sd)
                fit$ci.upper <- exp(outtab[5,  ] / Z.sd)
                fit$prob <- 1 - pchisq(outtab[6,  ], 1)
          } else {
                fit$method.ci <- "Wald"
                fit$ci.lower <- exp(coefs + qnorm(alpha/2) * vars^0.5)
                fit$ci.upper <- exp(coefs + qnorm(1 - alpha/2) * vars^0.5)
                fit$prob <- 1 - pchisq((coefs^2/vars), 1)
        }
        names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- cov.name
        attr(fit, "class") <- c("coxphf", "coxph")
        fit
}


"decomposeSurv" <-
function(
 formula, 
 data,
 sort=FALSE
)
### decomposes complex survival formula
### with time-factor interactions
### 2006-10
{
	orig.formula <- formula

	## expand formula if needed:
	repeat {
		terms <- terms(formula, data=data)
		fac <- attr(terms, "factors")
		needed <- rownames(fac)[!rownames(fac) %in% colnames(fac)][-1]

		if(length(needed) == 0) break

		formula <- as.formula(paste(as.character(formula)[2], "~", 
							as.character(formula)[3], "+",
							paste(needed, sep="+")))
	}

	## construct 3-col response:
	resp <- model.extract(model.frame(formula, data = data), response)
	if(ncol(resp) == 2)
		resp <- cbind(start=rep(0, nrow(resp)), resp)

	## sortieren nach STOPzeit und -Cens
	if(sort) {
	      sort <- order(resp[, 2],  -resp[, 3])
      	data <- data[sort, , drop=FALSE]
		resp <- resp[sort, ]
	}
		
	mm <- model.matrix(formula, data = data) ## Model-Matrix
	mm1 <- mm[, -1, drop=FALSE]	# w/o intercept
	
	terms <- terms(formula, data=data)
	fac <- attr(terms, "factors")
	labels <- attr(terms, "term.labels")
	
	## splittes by special chars
	f <- function(str)
		for(chars in c("(", ")", ":", " ", ",", "*", "^"))
			str <- unlist(strsplit(str, split=chars, fixed=TRUE))
		
	rowSplit <- sapply(rownames(fac), f)	# splitted effects 
	stopName <- tail(rowSplit[[1]], 2)[1]	# name of stoptime
	rowInter <- unlist(lapply(rowSplit[-1], function(z) any(z == stopName)))
		
	fac <- fac[-1, , drop=FALSE]	# omit Surv

	colSplit <- sapply(colnames(fac), f)
	colInter <- unlist(lapply(colSplit, function(z) any(z == stopName)))

	nTimes <- colSums(fac[rowInter, , drop=FALSE])
	nFac   <- colSums(fac[!rowInter, , drop=FALSE])
		
	inters <- (nFac>0) & (nTimes>0)
	NTDE <- sum(inters)

	
	timedata <- matrix(0, nrow(data), 0)
	timeind <- c()

	## loop for (time x effect)
	for(i in which(inters)) {
		## search pure time:
		ind <- (colSums(fac[rowInter, i] != fac[rowInter, , drop=FALSE]) == 0) & (nFac==0)
		timedata <- cbind(timedata, mm1[, ind, drop=FALSE])

		## search pure effect:
		ind <- (colSums(fac[!rowInter, i] != fac[!rowInter, , drop=FALSE]) == 0) & (nTimes == 0)
		timeind <- c(timeind, which(ind[!colInter]))
	}
	mm1 <- mm1[, !colInter, drop=FALSE]
	
	covnames <- c(colnames(mm1),
		paste(colnames(timedata), colnames(mm1)[timeind], sep=":")
	)

	## indicator to identify the original formula:
	ind <- covnames %in% colnames(attr(terms(orig.formula, data=data), "factors"))
	
	list(NTDE=NTDE, 			# number time dep. effects
		fac=fac, 			# factor matrix ..
		resp=resp, 			# N x 3 - response matrix
		mm1=mm1, 			# model matrix without time effects
		timedata=timedata, 	# matrix with time functions as columns
		timeind=timeind, 		# indicator of time-dependend effect
		covnames=covnames,	# names of covariates
		ind=ind			# indicator if some terms of not in formula
	)
}

"coxphfplot" <-
function(
 formula=attr(data, "formula"), 
 data=sys.parent(), 
 profile,      # righthand formula des zu plottenden Term (z.B. ~B oder ~A:D)
 pitch=.05,    # distances between points in std's
 limits,       # vector of MIN & MAX in std's, default=extremes of both CI's +-0.5 std. of beta
 alpha=.05, 
 maxit = 50, 
 maxhs = 5, 
 epsilon = 1e-006, 
 maxstep = 2.5, 
 firth = TRUE, 
 legend="center",	# see R-help of function <legend>, "" produces no legend
 ...			# other parameters to <legend>
) {
### by MP und GH, 2006-10
###
	## aufruf coxphf:
	fit <- coxphf(formula=formula, data=data, alpha=alpha, maxit=maxit, maxhs=maxhs, 
		epsilon=epsilon, maxstep=maxstep, firth=firth, pl=TRUE)
	coefs <- coef(fit)           ## "normale" Koeffizienten
	covs <- fit$var              ## Varianzen
	n <- nrow(data)
	
	obj <- decomposeSurv(formula, data, sort=TRUE)
      NTDE <- obj$NTDE
      mmm <- cbind(obj$mm1, obj$timedata)
	
	cov.name <- obj$covnames
	
      k <- ncol(obj$mm1)          # number covariates
      ones <- matrix(1, n, k+NTDE)
          
      ## standardisierung
	sd1 <- sd(obj$mm1)
	sd2 <- sd(obj$timedata)
	Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
	obj$mm1 <- scale(obj$mm1, FALSE, sd1)
	obj$timedata <- scale(obj$timedata, FALSE, sd2)
	mmm <- cbind(obj$mm1, obj$timedata)

	CARDS <- cbind(obj$mm1, obj$resp, ones, obj$timedata)   
      PARMS <- c(n, k, firth, maxit, maxhs, maxstep, epsilon, 0, 0, 0, 0, 0, 0, NTDE)
 
	#--> nun Berechnungen fuer Schleife
	formula2 <- as.formula(paste(as.character(formula)[2], 
                                        as.character(profile)[2], sep="~"))
      obj2 <- decomposeSurv(formula2, data, sort = TRUE)
      cov.name2 <- obj2$covnames       # Labels der Test-Fakt.
      pos <- match(cov.name2, cov.name) ## Position der Testfakt.

	std.pos <- diag(fit$var)[pos]^.5
	if(missing(limits)) {
		lim.pl <- (c(log(fit$ci.lower[pos]), log(fit$ci.upper[pos])) - 
		coef(fit)[pos]) / std.pos
		limits <- c(min(qnorm(alpha / 2), lim.pl[1]) - .5,
		max(qnorm(1 - alpha / 2), lim.pl[2]) + .5)
	}
	limits <- c(floor(limits[1]/pitch)*pitch, ceiling(limits[2]/pitch)*pitch)
	knots <- seq(limits[1], limits[2], pitch)
	iflag <- rep(1, k+NTDE)
	iflag[pos] <- 0
	offset <- rep(0, k+NTDE)
	nn <- length(knots)
	res <- matrix(knots, nn, 3) #initialisiere Werte
	dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
	for(i in 1:nn) {
		res[i, 2] <- coefs[pos] + covs[pos,pos]^.5 * knots[i]
		offset[pos] <- res[i, 2] * Z.sd[pos]  
		IOARRAY <- rbind(iflag, offset, matrix(0, 1+k+NTDE, k + NTDE))
		
		# --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
		value <- .Fortran("firthcox",
			CARDS,
			outpar = PARMS,
			outtab = IOARRAY, PACKAGE=coxphf)
		if(value$outpar[8])
			warning("Error in routine FIRTHCOX; parms8 <> 0")
		res[i, 3] <- value$outpar[11]
	}

	#### Graphischer Output:
	my.par <- act.par <- par()
	my.par$mai[3] <- 1.65 * act.par$mai[3]
###if(legend != "") 
###	my.par$mai[1] <- 2.00 * act.par$mai[1]
	par(mai=my.par$mai)
	ind <- (1:nn)[round(4*res[,1])==round(4*res[,1], 10)]
	if (length(ind)==0) ind <- 1:nn
	pp <- max(res[,3]) - .5*res[,1]^2
	plot(res[,-1], type="l", xlab=paste("BETA of", cov.name2)) ##Profile likelihood
	##lines(res[,2], pp, lty=4)  #<<<Wald approximative profile lik. >>>
	points(res[res[,1]==0,2], max(res[,3])) ##Maximum of likelihood
	segments(min(res[,2]), max(res[,3])-.5*qchisq(1-alpha, 1), 
		max(res[,2]), max(res[,3])-.5*qchisq(1-alpha, 1), lty=3) ##refer.line
	yy <- par("usr")[4] - (par("usr")[4]-par("usr")[3])*c(.9, .95)
	segments(coef(fit)[pos] - qnorm(alpha/2) * std.pos, yy[1], 
		coef(fit)[pos] - qnorm(1 - alpha/2) * std.pos, yy[1], lty=6) ##Wald-CI
	segments(log(fit$ci.lower[pos]), yy[2], 
		log(fit$ci.upper[pos]), yy[2], lty=8) ## prof.pen.lik.-CI
	axis(side=3, at=res[ind,2], labels=res[ind,1])
	mtext("distance from maximum in standard deviations", side=3, line=3)
	if(legend != "")
		legend(legend,
			##x=par("usr")[1], 
			##y=par("usr")[3]-.24*(par("usr")[4] - par("usr")[3]), 
			legend=c("Profile penalized likelihood", 
				paste(100*(1-alpha), "% - reference line", sep=""), 
				"Maximum of Likelihood",
				"Wald - confidence interval", 
				"Profile penalized likelihood CI"),
			lty=c(1, 3, -1, 6, 8),
			pch=c(-1, -1, 1, -1, -1),
			bty="n", ...
		)
	par(mai=act.par$mai)
	title("Profile likelihood")
	invisible(res)
}

"coxphftest" <-
function(
 formula=attr(data, "formula"), 
 data=sys.parent(), 
 test=~.,    # righthand formula fuer zu testende Faktoren (z.B. ~ B+D )
 values,     # fixierte Werte fuer die Faktoren, default=0 (z.B. c(2,0))
 maxit=50, 
 maxhs=5, 
 epsilon=1e-6, 
 maxstep=2.5, 
 firth=TRUE
) {
### MP, GH, 2006-10
###
	n <- nrow(data)
	obj <- decomposeSurv(formula, data, sort = TRUE)
    	NTDE <- obj$NTDE
    	mmm <- cbind(obj$mm1, obj$timedata)
    	cov.name <- obj$covnames
    	k <- ncol(obj$mm1)
    	ones <- matrix(1, n, k + NTDE)
	
	## standardisierung
	sd1 <- sd(obj$mm1)
	sd2 <- sd(obj$timedata)
	Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
	obj$mm1 <- scale(obj$mm1, FALSE, sd1)
	obj$timedata <- scale(obj$timedata, FALSE, sd2)
	mmm <- cbind(obj$mm1, obj$timedata)
          
	CARDS <- cbind(obj$mm1, obj$resp, ones, obj$timedata)   
	PARMS <- c(n, k, firth, maxit, maxhs, maxstep, epsilon, 0, 0, 0, 0, 0, 0, NTDE)
      IOARRAY <- rbind(rep(1, k+NTDE), matrix(0, 2+k+NTDE, k + NTDE))
      if(NTDE>0)
      	IOARRAY[4,(k+1):(k+NTDE)] <- obj$timeind
	storage.mode(CARDS) <- "double"
	storage.mode(PARMS) <- "double"
	storage.mode(IOARRAY) <- "double" #
	# --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
	value <- .Fortran("firthcox",
		CARDS,
		outpar = PARMS,
		outtab = IOARRAY, PACKAGE=coxphf)
	 if(value$outpar[8])
                warning("Error in routine FIRTHCOX; parms8 <> 0")
      loglik <- c(NA, value$outpar[11])

	############## now run test formula ################
	formula2 <- as.formula(paste(as.character(formula)[2], 
					as.character(test)[2], sep="~"))
	obj2 <- decomposeSurv(formula2, data, sort = TRUE)

	cov.name2 <- obj2$covnames[obj2$ind]	 # Labels der Test-Fakt.
	k2 <- length(cov.name2)		# Anzahl Faktoren
	if(!all(cov.name2 %in% cov.name))
		stop("formula test is not a subset of whole formula!")
	pos <- match(cov.name2, cov.name) ## Position der Testfakt.
	IOARRAY[1, pos] <- 0
	if(!missing(values))
		IOARRAY[2, pos] <- values * Z.sd[pos]
	# --------------- Aufruf Fortran - Makro FIRTHCOX ----------------------------------
	value <- .Fortran("firthcox",
		CARDS,
		outpar = PARMS,
		outtab = IOARRAY, PACKAGE=coxphf)
	if(value$outpar[8])
                warning("Error in routine FIRTHCOX; parms8 <> 0")
	loglik[1] <- value$outpar[11]
	
	testcov <- IOARRAY[2, ]
	testcov[-pos] <- NA
	names(testcov) <- cov.name

	## return object of class "coxphftest"
	fit <- list(testcov = testcov, loglik = loglik, df = k2, 
		prob = 1 - pchisq(2 * diff(loglik), k2),call = match.call())
	if(firth)
		fit$method <- "Penalized ML"
	else fit$method <- "Standard ML"
	attr(fit, "class") <- "coxphftest"
	fit
}

"print.coxphf" <-
function(
  x,   	# object of class coxphf
  ...		# dummy
)
### MP and GH
### 2006-10
{
	print(x$call)
	cat("Model fitted by", x$method)
	cat("\nConfidence intervals and p-values by", x$method.ci, "\n\n")
	
	out <- cbind(x$coefficients, diag(x$var)^0.5, exp(x$coefficients), 
		x$ci.lower, x$ci.upper, qchisq(1 - x$prob, 1), x$prob)
	dimnames(out) <- list(names(x$coefficients), c("coef", "se(coef)", 
		"exp(coef)", paste(c("lower", "upper"), 1 - x$alpha), "z", "p"))	

	if (x$method.ci != "Wald") 
		dimnames(out)[[2]][6] <- "Chisq"
	print(out)
	
	LL <- 2 * diff(x$loglik)
	cat("\nLikelihood ratio test=", LL, " on ", x$df, 
		" df, p=", 1 - pchisq(LL, x$df), ", n=", x$n, "\n\n", sep = "")

	invisible(x)
}

"print.coxphftest" <-
function(
  x,		# object of class coxphftest
  ...		# dummy
)
### MP, GH, 2006-10
{
	print(x$call)
	cat("Model fitted by", x$method, "\n\nFactors fixed as follows:\n")

	## output of factors and fixing values:
	print(x$testcov)
	LL <- 2 * diff(x$loglik) 
	out <- c(x$loglik[1], x$loglik[2], LL / 2)
	names(out) <- c("Restricted model", "Full model", "difference")

	## output of likelihoods:
	cat("\nLikelihoods:\n")
	print(out)

	## result summary:
	cat("\nLikelihood ratio test=", LL, " on ", x$df, 
		" df, p=", x$prob, "\n", sep="")
	
	invisible(x)
}

"summary.coxphf" <-
function(
  object,		# object of class coxphf
  ...			# dummy
)
### MP and GH
### 2006-10
{
	print(object$call)
	cat("\nModel fitted by", object$method)
	cat("\nConfidence intervals and p-values by", object$method.ci, "\n\n")

	out <- cbind(object$coefficients, diag(object$var)^0.5, 
		exp(object$coefficients), 
		object$ci.lower, object$ci.upper,  
		qchisq(1 - object$prob, 1), object$prob)
	dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)", 
		"exp(coef)", paste(c("lower", "upper"), 1 - object$alpha), "z", "p"))
	if (object$method.ci != "Wald") 
		dimnames(out)[[2]][6] <- "Chisq"

	print(out)

	LL <- 2 * diff(object$loglik)
	cat("\nLikelihood ratio test=", LL, " on ", object$df, 
		" df, p=", 1 - pchisq(LL, object$df), ", n=", object$n, sep = "")
	wald.z <- t(coef(object)) %*% solve(object$var) %*% coef(object) 
	cat("\nWald test =", wald.z, "on", object$df, "df, p =", 
		1 - pchisq(wald.z, object$df))
	cat("\n\nCovariance-Matrix:\n")
	print(object$var)
	
	invisible(object)
}

