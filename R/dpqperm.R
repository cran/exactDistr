dperm <- function(x, scores, m)
{
	if (x < 0) return(0)
	eq <- equiscores(scores)
	cp <- cperm(eq, m)
	
	# why does RVAL <- cp$Prob[cp$T == x] not work ???
	
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
	if (length(q) != 1) stop("x is not a real number")
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
	if (!is.integer(scores)) 
	{
		fscore <- scores - floor(scores)
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

	# Shift algorithmus for the exact Wilcoxon distribution
	# by B. Streitberg and J. Roehmel: Exact Distributions for
	# Permutations and Rank Tests: An Introduction to Some Recently
	# Published Algorithms, Statistical Software Newsletter 1984, 
	# Vol. 12, No. 1, pp. 10-17
	#   
	# Is equal to dwilcox if scores are untied. 
	# Can be extended to other rank tests by modificated scores

	a <- escores$scores 

	N <- length(a)
	H <- matrix(1)
	c <- sum(a[(m+1):N])

	for (i in 1:N)
	{
		if (a[i] > 0)
		{
			A <- cbind(H, matrix(rep(0, a[i]*nrow(H)), ncol=a[i]))
			A <- rbind(A, 0)
			B <- cbind(matrix(rep(0, a[i]*nrow(H)), ncol=a[i]), H)
			B <- rbind(0, B)
			H <- A + B
			if (nrow(H) > m+1)
				H <- H[1:(m+1), ]
			if (ncol(H) > c+1)
				H <- H[,1:(c+1)]
		} else {
			H <- H + H
		}			
	}
	prob <- H[m+1,]
	t <- which(prob != 0)
	prob <- prob[t]
	fact <- gamma(m+1)*gamma(N -m +1)
	fact <- fact/gamma(N+1)
	prob <- fact*prob
	t <- t - 1
	t <- 1/escores$factor*t

	# return the density P( T = t) = prob for all possible t

	RVAL <- list(T = t, Prob = prob)
	class(RVAL) <- "cperm"
	return(RVAL)
}

	