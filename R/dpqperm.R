
dperm <- function(x, scores, m)
{
	if (x < 0) return(0)
	eq <- equiscores(scores)
	cp <- cperm(eq, m)
	
	RVAL <- cp$Prob[cp$T == x[1]]
	if (length(RVAL) == 0) RVAL <- 0
	if (length(x) > 1)
	{
		for (i in c(2:length(x)))
		{
			pr <- cp$Prob[cp$T == x[i]]
			if (length(pr) == 0) pr <- 0
			RVAL <- c( RVAL, pr)
		}
	} 
	if (is.null(RVAL)) RVAL <- rep(0, length(x))
	return(RVAL) 
}

pperm <- function(q, scores, m)
{
	if (length(q) != 1) stop("q is not a real number")
        if (q < 0) return(0)
        eq <- equiscores(scores)
        cp <- cperm(eq, m)
	return(sum(cp$Prob[cp$T <= q]))
}

qperm <- function(p, scores, m)
{
	if (p < 0 || p > 1) stop("p is not a probability")
	eq <- equiscores(scores)                          
        cp <- cperm(eq, m)
	return(max(cp$T[cumsum(cp$Prob) < p]) + 1)
}

equiscores <- function(scores)
{
	fscore <- scores - floor(scores)
	
	if (!all(fscore == 0)) 
	{
		if (all(fscore[fscore != 0] == 0.5))
		{
			factor <- 2
			scores <- factor*scores
		} else {
			stop("midrank scores allowed only!")
		}
	} else { 
		factor <- 1
	}

	RVAL <- list(scores = scores, factor = factor)
	class(RVAL) <- "equis"
	return(RVAL)
}


cperm <- function(escores, m)
{

	if (!(class(escores) == "equis"))
		stop("scores are not of class equis") 

	N <- length(escores$scores)

	prob <- rep(0, max(cumsum(escores$scores)))

	if (N == m) {

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

	t <- 1/escores$factor*t

	RVAL <- list(T = t, Prob = prob)
	class(RVAL) <- "cperm"
	return(RVAL)
}

	