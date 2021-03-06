
toltest <- function(x, scores, m)
    sum(rev(sort(abs(scores*x - round(scores*x))))[1:m])/x

toltest <- function(x, scores, m) {
    # something more sensitive
    a <- sort(scores*x - round(scores*x))
    upper <- a[1:m]       
    upper <- sum(abs(upper[upper < 0]))
    lower <- a[(length(a) - m):length(a)]
    lower <- sum(abs(lower[lower > 0]))
    max(upper, lower)/x   
}


dperm <- function(x, scores, m, paired = NULL, tol = 0.01, fact = NULL)
{
    eq <- equiscores(scores, m, tol, fact)
    cp <- cperm(eq, m, paired)
    RVAL <- rep(0, length(x))
    RVAL[x %in% cp$T ] <- cp$Prob[cp$T %in% x]
    return(RVAL) 
}

pperm <- function(q, scores, m, paired = NULL, tol = 0.01, fact = NULL)
{
    eq <- equiscores(scores, m, tol, fact)
    cp <- cperm(eq, m, paired)
    RVAL <- c()
    for (i in q) {
        prob <- cp$Prob[cp$T <= i]
        ifelse(length(prob) >= 1, RVAL <- c(RVAL, sum(prob)),
                                  RVAL <- c(RVAL, 0))
    }
    return(RVAL)
}

qperm <- function(p, scores, m, paired = NULL, tol = 0.01, fact = NULL)
{
    if (any(p < 0) || any(p > 1)) {
        warning("p is not a probability")
        return(NaN)
    }
    eq <- equiscores(scores, m, tol, fact)                          
    cp <- cperm(eq, m, paired)
    cs <- cumsum(cp$Prob)
    RVAL <- c()
    for (i in p) {
        quant <- which(cs < i)
        if (length(quant) == 0) quant <- 0
        RVAL <- c(RVAL, cp$T[max(quant) + 1])
    }
    return(RVAL)
}

rperm <- function(size, scores, m)
    sapply(1:size, dummy <- function(x) sum(sample(scores,m)))

equiscores <- function(scores, m, tol = 0.01, fact=NULL)
{
    fscore <- scores - floor(scores)
    
    if (all(fscore == 0))
    { 
        # integer valued scores
        fact <- 1
        add <- min(scores) - 1
	scores <- scores - add
    } else {
        if (all(fscore[fscore != 0] == 0.5))
        {
            # midranked scores
            fact <- 2
            scores <- scores*fact
            add <- min(scores) - 1
            scores <- scores - add
        } else {
            # rational scores
            ssc <- sort(scores)
            b <- ceiling(1/min(ssc[2:length(ssc)] - ssc[1:(length(ssc)-1)]))
            if (is.null(fact) || fact < b ) {
                if (toltest(b, scores, m) <= tol)
                    fact <- b
                else {
                    # do not induce more than 20.000 columns
                    maxfact <- min(1000, round(20000/sum(scores -
                                                    min(scores))))   
                    if (maxfact < b)
                        fact <- maxfact
                    else {
                        test <- function(x, tol, sc)
                                ifelse(toltest(x, sc, m) - tol > 0, 1/x, x)
                        fact <- optim(10, test, tol=tol, sc=scores)$par
                        if (fact > maxfact) {
                            fact <- maxfact
                            warning(paste("cannot hold tol, tolerance:",
                                     round(toltest(maxfact, scores, m), 6)))
                        } else fact <- min(fact)
                    }
                }
            } 
            scores <- round(scores * fact)
            add <- min(scores)-1
            scores <- scores - add
        }
    }

    RVAL <- list(scores = scores, fact = fact, add = add)
    class(RVAL) <- "equis"
    return(RVAL)
}


cperm <- function(escores, m, paired = NULL)
{
    if (!(class(escores) == "equis"))
        stop("scores are not of class equis") 

    N <- length(escores$scores)

    prob <- rep(0, max(cumsum(escores$scores)))

    if (is.null(paired))
        paired <- (N == m)
    else 
        paired <- (N == m) && paired  

    if (paired) {
        # paired two sample situation
        prob <- c(0, prob)
        prob <- .C("cpermdist1", prob = as.double(prob),
                   as.integer(escores$scores), as.integer(N))$prob
        t <- which(prob != 0)
        prob <- prob[t]
        # 0 is possible
	t <- t - 1
    } else {
        # independent samples
        col <- sum(sort(escores$scores)[(N + 1 - m):N])
        scores <- rep(1, N)
        prob <- .C("cpermdist2", prob = as.double(prob), as.integer(m),
                as.integer(col), as.integer(scores), as.integer(escores$scores),
                as.integer(N), as.integer(1))$prob
        t <- which(prob != 0)
        prob <- prob[t]
    }
    t <- (t + escores$add*m)/escores$fact
    RVAL <- list(T = t, Prob = prob)
    class(RVAL) <- "cperm"
    return(RVAL)
}
    