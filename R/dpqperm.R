
dperm <- function(x, scores, m, paired = NULL, fact = NULL)
{
    if (x < 0) return(0)
    eq <- equiscores(scores, fact)
    cp <- cperm(eq, m, paired)
    RVAL <- rep(0, length(x))
    RVAL[x %in% cp$T ] <- cp$Prob[cp$T %in% x]
    return(RVAL) 
}

pperm <- function(q, scores, m, paired = NULL, fact = NULL)
{
    if (length(q) != 1) stop("q is not a real number")
    if (q < 0) return(0)
    eq <- equiscores(scores, fact)
    cp <- cperm(eq, m, paired)
    return(sum(cp$Prob[cp$T <= q]))
}

qperm <- function(p, scores, m, paired = NULL, fact = NULL)
{
    if (p < 0 || p > 1) stop("p is not a probability")
    eq <- equiscores(scores, fact)                          
    cp <- cperm(eq, m)
    cp$T[max(which(cumsum(cp$Prob) < p)) + 1]
}

equiscores <- function(scores, fact=NULL)
{
    fscore <- scores - floor(scores)
    
    if (all(fscore == 0))
    { 
        # integer valued scores
        fact <- 1
        add <- min(scores) - 1
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
            b <- 2/min(ssc[2:length(ssc)] - ssc[1:(length(ssc)-1)])
            if (is.null(fact) || fact < b )
            fact <- b
            scores <- round(scores * fact)
            add <- min(0, min(scores)-1)
            scores <- scores - add
        }
    }

    RVAL <- list(scores = scores, fact = fact, add = add)
    class(RVAL) <- "equis"
    return(RVAL)
}


cperm <- function(escores, m, paired = NULL, fact = NULL)
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
        prob <- .C("cpermdist1", prob = as.double(prob), as.integer(escores$scores), as.integer(N))$prob
        t <- which(prob != 0)
        prob <- prob[t]
        # 0 is possible
	t <- t - 1
    } else {
        # independent samples
        col <- sum(sort(escores$scores)[(N + 1 - m):N])
        scores <- rep(1, N)
        prob <- .C("cpermdist2", prob = as.double(prob), as.integer(m),as.integer(col), as.integer(scores), as.integer(escores$scores), as.integer(N))$prob
        t <- which(prob != 0)
        prob <- prob[t]
    }
    t <- (t + escores$add*m)/escores$fact
    RVAL <- list(T = t, Prob = prob)
    class(RVAL) <- "cperm"
    return(RVAL)
}
    