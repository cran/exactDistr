\name{dperm}
\alias{dperm}
\title{Distribution of Permutation Tests}
\description{
Density, distribution function and quantile function 
for the distribution of Permutation Tests using the Shift-Algorithm by
Streitberg and R\"ohmel.
}
\usage{
dperm(x, scores, m, paired=NULL, fact=NULL)
pperm(q, scores, m, paired=NULL, fact=NULL)
qperm(p, scores, m, paired=NULL, fact=NULL)
}
\arguments{
  \item{x, q}{ vector of quantiles }
  \item{p}{vector of probabilities}
  \item{scores}{the (mid)ranks of the observations of the \code{x} (first
\code{m} elements) and \code{y} sample. }
  \item{m}{ sample size of the \code{x} sample. If \code{m = length(x)}
scores of paired observations are assumed.}
  \item{paired}{ logical. Indicates if paired observations are used. Only
needed for the distribution of the total sum of the scores (which has mass 1
at the point \code{sum(scores)}}
  \item{fact}{ rational scores are mapped into integers by multiplication.
Use \code{fact} instead of an internal factor.}
}
\details{
  The exact distribution of the sum of the first \code{m} scores is
evaluated using the Shift-Algorithm by Streitberg and R\"ohmel under the
hypothesis of exchangeability.
}
\value{
\code{dperm} gives the density, \code{pperm} gives the distribution
     function and \code{qperm} gives the quantile function.

}
\references{ 

	Bernd Streitberg and Joachim R\"ohmel (1986). 
	Exact Distributions For Permutations and Rank Tests: 
	An Introduction to Some Recently Published Algorithms 
	Statistical Software Letter, 12 No. 1, 10-17.

	Bernd Streitberg and Joachim R\"ohmel (1987)
	Exakte Verteilungen f\"ur Rang- und Randomisierungstests 
	im allgemeinen $c$-Stichprobenfall
	EDV in Medizin und Biologie, 18 No. 1, 12-19 
}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

# exact p-value of the Wilcoxon test for a tied sample

x <- c(0.5, 0.5, 0.6, 0.6, 0.7, 0.8, 0.9)
y <- c(0.5, 1.0, 1.2, 1.2, 1.4, 1.5, 1.9, 2.0)
r <- rank(c(x,y))
pexact <- pperm(sum(r[seq(along=x)]), r, 7)

# Compare the exact algorithm as implemented in ctest and the
# Streitberg-Roehmel for untied samples
 
# Wilcoxon:

n <- 10
x <- rnorm(n, 2)
y <- rnorm(n, 3)
r <- rank(c(x,y))

# exact distribution using Streitberg-Roehmel

dwexac <- dperm((n*(n+1)/2):(n^2 + n*(n+1)/2), r, n)
su <- sum(dwexac)                             # should be something near 1 :-)
if (su != 1) stop("dum(dwexac) not equal 1")

# exact distribution using dwilcox

dw <- dwilcox(0:(n^2), n, n)

# compare the two distributions:

plot(dw, dwexac)      # should give a "perfect" line

# Wilcoxon signed rank test

x <- rnorm(10, 5)
y <- rnorm(10, 5)
r <- rank(abs(x - y))
pperm(sum(r[x - y > 0]), r, length(r))
wilcox.test(x,y, paired=T, alternative="less")
psignrank(sum(r[x - y > 0]), length(r))

# Ansari-Bradley

n <- 10
x <- rnorm(n, 2, 1)
y <- rnorm(n, 2, 2)

# exact distribution using Streitberg-Roehmel

r <- rank(c(x,y))
sc <- pmin(r, 2*n - r +1)
dabexac <- dperm(0:(n*(2*n+1)/2), sc, n)
sum(dabexac)
tr <- which(dabexac > 0)

# exact distribution using dansari (wrapper to ansari.c in ctest)

dab <- dansari(0:(n*(2*n+1)/2), n, n)

# compare the two distributions:

plot(dab[tr], dabexac[tr])

# rational scores are allowed (but only result in an approximation)
# e.g. v.d. Waeren test

x <- rnorm(10)
y <- rnorm(10)
N <- length(x) + length(y)
r <- rank(c(x,y))
scores <- qnorm(r/(N+1))
X <- sum(scores[seq(along=x)])  # <- v.d. Waeren normal quantile statistic

# critical value, two-sided test
abs(qperm(0.025, scores, length(x), fact=16))

# p-values

p <- pperm(abs(X), scores, length(x))
min(p, 1-p)

# check if [qp]wilcox and [pq]ansari are equal with [qp]perm

source(system.file("test1.R", pkg="exactDistr"))
source(system.file("test2.R", pkg="exactDistr"))

}
\keyword{exact distribution}
