wilcox.exact <-
function(x, y = NULL, alternative = c("two.sided", "less", "greater"), 
         mu = 0, paired = FALSE, exact = NULL, correct = TRUE,
         conf.int = FALSE, conf.level = 0.95) 
{
    alternative <- match.arg(alternative)
    if(!missing(mu) && ((length(mu) > 1) || !is.finite(mu))) 
        stop("mu must be a single number")
    if(conf.int) {
        if(!((length(conf.level) == 1)
             && is.finite(conf.level)
             && (conf.level > 0)
             && (conf.level < 1)))
            stop("conf.level must be a single number between 0 and 1")
    }

    if(!is.null(y)) {
        DNAME <- paste(deparse(substitute(x)), "and",
                       deparse(substitute(y)))
        if(paired) {
            if(length(x) != length(y)) 
                stop("x and y must have the same length")
            OK <- complete.cases(x, y)
            x <- x[OK] - y[OK]
            y <- NULL
        }
        else {
            x <- x[is.finite(x)]
            y <- y[is.finite(y)]
        }
    } else {
        DNAME <- deparse(substitute(x))
        if(paired) 
            stop("y missing for paired test")
        x <- x[is.finite(x)]
    }

    if(length(x) < 1) 
        stop("not enough (finite) x observations")
    CORRECTION <- 0
    if(is.null(y)) {
        METHOD <- "Wilcoxon signed rank test"
        x <- x - mu
        ZEROES <- any(x == 0)
        if(ZEROES) 
            x <- x[x != 0]
        n <- length(x)
        if(is.null(exact)) 
            exact <- (n < 50)
        r <- rank(abs(x))
        STATISTIC <- sum(r[x > 0])
        names(STATISTIC) <- "V"
        TIES <- (length(r) != length(unique(r)))

	### TEST ### psignrank ist gleich schnell wie pperm

        if (FALSE) { # if(exact && !TIES && !ZEROES) {
            PVAL <-
                switch(alternative,
                       "two.sided" = {
                           if(STATISTIC > (n * (n + 1) / 4)) 
                               p <- 1 - psignrank(STATISTIC - 1, n)
                           else p <- psignrank(STATISTIC, n)
                           min(2 * p, 1)
                       },
                       "greater" = 1 - psignrank(STATISTIC - 1, n),
                       "less" = psignrank(STATISTIC, n))
            if(conf.int && !is.na(x)) {
                ## Exact confidence intervale for the median in the
                ## one-sample case.  When used with paired values this
                ## gives a confidence interval for mean(x) - mean(y).
                x <- x + mu             # we want a conf.int for the median
                alpha <- 1 - conf.level
                diffs <- outer(x, x, "+")
                diffs <- sort(diffs[!lower.tri(diffs)]) / 2
                cint <-
                    switch(alternative,
                           "two.sided" = {
                               qu <- qsignrank(alpha / 2, n)
                               if(qu == 0) qu <- 1
                               ql <- n*(n+1)/2 - qu
                               uci <- diffs[qu]
                               lci <- diffs[ql+1]
                               c(uci, lci)        
                           },
                           "greater"= {
                               qu <- qsignrank(alpha, n)
                               if(qu == 0) qu <- 1
                               uci <- diffs[qu]
                               c(uci, NA)
                           },
                           "less"= {
                               qu <- qsignrank(alpha, n)
                               if(qu == 0) qu <- 1
                               ql <- n*(n+1)/2 - qu
                               lci <- diffs[ql+1]
                               c(NA, lci)        
                           })
                attr(cint, "conf.level") <- conf.level    
            }
        } else {
            if (exact) {
                PVAL <-
                   switch(alternative,
                          "two.sided" = {
                              if(STATISTIC > (n * (n + 1) / 4)) 
                                  p <- 1 - pperm(STATISTIC - 1, r, n)
                              else p <- pperm(STATISTIC, r, n)
                              min(2 * p, 1)
                          },
                          "greater" = 1 - pperm(STATISTIC - 1, r, n),
                          "less" = pperm(STATISTIC, r, n))
                if(conf.int && !is.na(x)) {
                    ## Exact confidence intervale for the median in the
                    ## one-sample case.  When used with paired values this
                    ## gives a confidence interval for mean(x) - mean(y).
                    x <- x + mu             # we want a conf.int for the median
                    alpha <- 1 - conf.level
                    diffs <- outer(x, x, "+")
                    diffs <- sort(diffs[!lower.tri(diffs)]) / 2
                    w <- cumsum(rep(1, max(cumsum(r))))
                    dup <- which(duplicated(diffs))
                    indx <- 1:length(diffs)
                    indx[dup] <- NA
                    w <- w[!is.na(indx)]
                    diffs <- diffs[!is.na(indx)]
                    cint <-
                        switch(alternative,
                               "two.sided" = {
                                   qu <- qperm(alpha/2, r, n) 
                                   ql <- qperm(1-alpha/2, r, n)
                                   if (qu <= min(w)) uci <- min(diffs)
                                   else uci <- max(diffs[w <= qu])
                                   if (ql >= max(w)) lci <- max(diffs)
                                   else lci <- min(diffs[w > ql])
                                   c(uci, lci)
                               },
                               "greater"= {
                                   qu <- qperm(alpha, r, n) 
                                   if (qu <= min(w)) uci <- min(diffs)
                                   else uci <- max(diffs[w <= qu])
                                   c(uci, NA)
                               },
                               "less"= {
                                   ql <- qperm(1-alpha, r, n)
                                   if (ql >= max(w)) lci <- max(diffs)
                                   else lci <- min(diffs[w > ql])
                                   c(NA, lci)
                         })
                     attr(cint, "conf.level") <- conf.level    
                }
            } else {
                NTIES <- table(r)
                z <- STATISTIC - n * (n + 1)/4
                SIGMA <- sqrt(n * (n + 1) * (2 * n + 1) / 24
                              - sum(NTIES^3 - NTIES) / 48)
                if(correct) {
                    CORRECTION <-
                        switch(alternative,
                               "two.sided" = sign(z) * 0.5,
                               "greater" = 0.5,
                               "less" = -0.5)
                    METHOD <- paste(METHOD, "with continuity correction")
                }

                PVAL <- pnorm((z - CORRECTION) / SIGMA)
                if(alternative == "two.sided") 
                    PVAL <- 2 * min(PVAL, 1 - PVAL)
                if(alternative == "greater") 
                    PVAL <- 1 - PVAL

                if(conf.int && !is.na(x)) {
                    ## Asymptotic confidence intervale for the median in the
                    ## one-sample case.  When used with paired values this
                    ## gives a confidence interval for mean(x) - mean(y).
                    ## Algorithm not published, thus better documented here.
                    x <- x + mu
                    alpha <- 1 - conf.level
                    ## These are sample based limits for the median
                    mumin <- min(x)
                    mumax <- max(x)
                    ## wdiff(d, zq) returns the abolute difference between
                    ## the asymptotic Wilcoxon statistic of x - mu - d and
                    ## the quantile zq 
                    wdiff <- function(d, zq) {
                        xd <- x  - d
                        xd <- xd[xd != 0]
                        nx <- length(xd)
                        dr <- rank(abs(xd))
                        zd <- sum(dr[xd > 0])
                        NTIES.CI <- table(dr)
                        zd <- zd - nx * (nx + 1)/4
                        SIGMA.CI <- sqrt(nx * (nx + 1) * (2 * nx + 1) / 24
                                         - sum(NTIES.CI^3 -  NTIES.CI) / 48)
                        if(correct) {
                            CORRECTION.CI <-
                                switch(alternative,
                                       "two.sided" = sign(z) * 0.5,
                                       "greater" = 0.5,
                                       "less" = -0.5)
                        }
                        zd <- (zd - CORRECTION.CI) / SIGMA.CI
                        abs(zd - zq)
                    }
                    ## Here we optimize the function wdiff in d over the set
                    ## c(mumin, mumax).
                    ##
                    ## This returns a value from c(mumin, mumax) for which
                    ## the asymptotic Wilcoxon statistic is equal to the
                    ## quantile zq.  This means that the statistic is not
                    ## within the critical region, and that implies that d
                    ## is a confidence limit for the median.
                    ##
                    ## As in the exact case, interchange quantiles.
                    cint <- switch(alternative, "two.sided" = {
                        u <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                      zq=qnorm(1-alpha/2))$minimum
                        l <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                      zq=qnorm(alpha/2))$minimum
                        c(u, l)
                    }, "greater"= {
                        u <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                      zq=qnorm(1-alpha))$minimum
                        c(u, NA)
                    }, "less"= {
                        l <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                      zq=qnorm(alpha))$minimum
                        c(NA, l)
                    })
                    attr(cint, "conf.level") <- conf.level    
                }

#                if(exact && TIES) {
#                warning("Cannot compute exact p-value with ties")
#                if(conf.int)
#                    warning(paste("Cannot compute exact confidence",
#                                  "interval with ties"))
#                }
#                if(exact && ZEROES) {
#                   warning("Cannot compute exact p-value with zeroes")
#                   if(conf.int)
#                       warning(paste("Cannot compute exact confidence",
#                                     "interval with zeroes"))
#                }
            }            
	}
    }
    else {
        if(length(y) < 1) 
            stop("not enough y observations")
        METHOD <- "Wilcoxon rank sum test"
        r <- rank(c(x - mu, y))
        n.x <- length(x)
        n.y <- length(y)
        if(is.null(exact)) 
            exact <- (n.x < 50) && (n.y < 50)
        STATISTIC <- sum(r[seq(along = x)]) - n.x * (n.x + 1) / 2
        names(STATISTIC) <- "W"
        TIES <- (length(r) != length(unique(r)))
        if(exact) {

	    ### TEST ### pperm ist langsamer als pwilcox

	    if (FALSE) { # if (!TIES) {
                PVAL <-
                    switch(alternative,
                           "two.sided" = {
                               if(STATISTIC > (n.x * n.y / 2)) 
                                   p <- 1 - pwilcox(STATISTIC - 1, n.x, n.y)
                               else
                                   p <- pwilcox(STATISTIC, n.x, n.y)
                               min(2 * p, 1)
                           },
                           "greater" = 1 - pwilcox(STATISTIC - 1, n.x, n.y), 
                           "less" = pwilcox(STATISTIC, n.x, n.y))
                if(conf.int) {
                    ## Exact confidence interval for the location parameter 
                    ## mean(y) - mean(x) in the two-sample case (cf. the
                    ## one-sample case).
                    alpha <- 1 - conf.level
                    diffs <- sort(outer(y, x, "-"))
                    cint <-
                        switch(alternative,
                               "two.sided" = {
                                   qu <- qwilcox(alpha/2, n.x, n.y)
                                   if(qu == 0) qu <- 1
                                   ql <- n.x*n.y - qu
                                   uci <- diffs[qu]
                                   lci <- diffs[ql + 1]
                                   c(uci, lci)
                               },
                               "greater"= {
                                   qu <- qwilcox(alpha, n.x, n.y)
                                   if(qu == 0) qu <- 1
                                   uci <- diffs[qu]
                                   c(uci, NA)
                               },
                               "less"= {
                                   qu <- qwilcox(alpha, n.x, n.y)
                                   if(qu == 0 ) qu <- 1
                                   ql <- n.x*n.y - qu
                                   lci <- diffs[ql + 1]
                                   c(NA, lci)
                               })
                    attr(cint, "conf.level") <- conf.level    
                }
            } else {
		# now the exact case using Streitberg/Roehmel
                PVAL <-
                    switch(alternative,
                           "two.sided" = {
                               if(STATISTIC > (n.x * n.y / 2)) 
                                   p <- 1 - pperm(STATISTIC - 1 + n.x*(n.x +1)/2, r, n.x)
                               else
                                   p <- pperm(STATISTIC + n.x*(n.x +1)/2, r, n.x)
                               min(2 * p, 1)
                           },
                           "greater" = 1 - pperm(STATISTIC - 1 + n.x*(n.x+1)/2, r, n.x), 
                           "less" = pperm(STATISTIC+ n.x*(n.x +1)/2, r, n.x))
                if(conf.int) {
                    ## Exact confidence interval for the location parameter 
                    ## mean(y) - mean(x) in the two-sample case (cf. the
                    ## one-sample case).
                    alpha <- 1 - conf.level
                    diffs <- sort(outer(y, x, "-"))
                    w <- cumsum(rep(1, n.x*n.y))
                    dup <- which(duplicated(diffs))
                    indx <- 1:length(diffs)
                    indx[dup] <- NA
                    w <- w[!is.na(indx)]
                    diffs <- diffs[!is.na(indx)]   
                    cint <-
                        switch(alternative,
                               "two.sided" = {
                                   qu <- qperm(alpha/2, r, n.x) - n.x*(n.x+1)/2
                                   ql <- qperm(1-alpha/2, r, n.x) - n.x*(n.x+1)/2
                                   if (qu <= min(w)) uci <- min(diffs)
                                   else uci <- max(diffs[w <= qu])
                                   if (ql >= max(w)) lci <- max(diffs)
                                   else lci <- min(diffs[w > ql])
                                   c(uci, lci)
                               },
                               "greater"= {
	                           qu <- qperm(alpha, r, n.x) - n.x*(n.x+1)/2
                                   if (qu <= min(w)) uci <- min(diffs)
                                   else uci <- max(diffs[w <= qu])
                                   c(uci, NA)
                               },
                               "less"= {
                                   ql <- qperm(1-alpha, r, n.x) - n.x*(n.x+1)/2
                                   if (ql >= max(w)) lci <- max(diffs)
                                   else lci <- min(diffs[w > ql]) 
                                   c(NA, lci)
                              })
                     attr(cint, "conf.level") <- conf.level    
                }
            }
        } else {
            NTIES <- table(r)
            z <- STATISTIC - n.x * n.y / 2
            SIGMA <- sqrt((n.x * n.y / 12) *
                          ((n.x + n.y + 1)
                           - sum(NTIES^3 - NTIES)
                           / ((n.x + n.y) * (n.x + n.y - 1))))
            if(correct) {
                CORRECTION <- switch(alternative,
                                     "two.sided" = sign(z) * 0.5,
                                     "greater" = 0.5,
                                     "less" = -0.5)
                METHOD <- paste(METHOD, "with continuity correction")
            }
            PVAL <- pnorm((z - CORRECTION)/SIGMA)
            if(alternative == "two.sided") 
                PVAL <- 2 * min(PVAL, 1 - PVAL)
            if(alternative == "greater") 
                PVAL <- 1 - PVAL

            if(conf.int) {
                ## Asymptotic confidence interval for the location
                ## parameter mean(y) - mean(x) in the two-sample case
                ## (cf. one-sample case).
                ##
                ## Algorithm not published, for a documentation see the
                ## one sample case.
                alpha <- 1 - conf.level
                mumin <- min(y) - max(x)
                mumax <- max(y) - min(x)
                wdiff <- function(d, zq) {
                    dr <- rank(c(x - mu, y - d)) 
                    NTIES.CI <- table(dr)
                    dz <- (sum(dr[seq(along = x)])
                           - n.x * (n.x + 1) / 2 - n.x * n.y / 2)
                    if(correct) {
                        CORRECTION.CI <-
                            switch(alternative,
                                   "two.sided" = sign(dz) * 0.5,
                                   "greater" = 0.5,
                                   "less" = -0.5)        
                    }
                    SIGMA.CI <- sqrt((n.x * n.y / 12) *
                                     ((n.x + n.y + 1)
                                      - sum(NTIES.CI^3 - NTIES.CI)
                                      / ((n.x + n.y) * (n.x + n.y - 1))))
                    dz <- (dz - CORRECTION.CI) / SIGMA.CI
                    abs(dz - zq)
                }
                cint <- switch(alternative, "two.sided" = {
                    u <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(alpha/2))$minimum
                    l <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(1 - alpha/2))$minimum
                    c(u, l)
                }, "greater"= {
                    u <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(alpha))$minimum
                    c(u, NA)
                }, "less"= {
                    l <- optimize(wdiff, c(mumin, mumax), tol=1e-4,
                                  zq=qnorm(1 - alpha))$minimum
                    c(NA, l)
                })
                attr(cint, "conf.level") <- conf.level    
            }

#            if(exact && TIES) {
#                warning("Cannot compute exact p-value with ties")
#                if(conf.int)
#                    warning(paste("Cannot compute exact confidence",
#                                  "intervals with ties"))
#            }
        }
    }

    RVAL <- list(statistic = STATISTIC,
                 parameter = NULL,
                 p.value = PVAL, 
                 null.value = c(mu = mu),
                 alternative = alternative,
                 method = METHOD, 
                 data.name = DNAME)
    if(conf.int)
        RVAL$conf.int <- cint
    class(RVAL) <- "htest"
    return(RVAL)
}
