
# independent observations

library(exactDistr)

hansi <- c()
seppl <- c()
for (i in 1:10)
{
	m <- sample(10:100, 1)
	if (runif(1) < 0.5)
		n <- sample(10:100, 1)
	else	
		n <- m
	score <- sample(n+m)
	val <- sample(0:(m*n), 1)
	# cat("m: ", m, "n: ", n, " val: ", val, "\n")
	hansi <- c(hansi,  pwilcox(val, m, n))
	cat("pwilcox: ", hansi[length(hansi)])
	seppl <- c(seppl, pperm(val + m*(m+1)/2, score, m))
	cat(" pperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")

hansi <- c()
seppl <- c()
for (i in 1:10)
{
        m <- sample(10:100, 1)
        if (runif(1) < 0.5)
                n <- sample(10:100, 1)
        else
                n <- m
        score <- sample(n+m)
        prob <- runif(1)
        # cat("m: ", m, "n: ", n, " prob: ", prob, "\n")
        hansi <- c(hansi,  qwilcox(prob, m, n))
        cat("qwilcox: ", hansi[length(hansi)])
        seppl <- c(seppl, qperm(prob, score, m) - m*(m+1)/2)
        cat(" qperm: ", seppl[length(seppl)], "\n")
}

cat("Max difference: ", max(abs(hansi - seppl)), "\n")
