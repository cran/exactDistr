\name{qperm}
\alias{qperm}
\title{Distribution of Permutation Tests}
\description{
Density, distribution function and quantile function 
for the distribution of the Permutation Tests using the Shift-Algorithm by
Streitberg and R\"ohmel.
}
\usage{
dperm(x, scores, m)
pperm(q, scores, m)
qperm(p, scores, m)
}
\arguments{
  \item{x, q}{ vector of quantiles }
  \item{p}{vector of probabilities}
  \item{scores}{the (mid)ranks of the observations of the \code{x} (first
\code{m} elements) and \code{y} sample. }
  \item{m}{ sample size of the \code{x} sample }
}
\details{
  The exact distribution of the sum of the first \code{m} scores is
evaluated using the Shift-Algorithm by Streitberg and R\"ohmel.
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

x <- rnorm(10, 2)
y <- rnorm(10, 3)
r <- rank(c(x,y))

# exact distribution using Streitberg-Roehmel

dwexac <- dperm(0:100, r, 10)
sum(dwexac)                             # should be something near 1 :-)
tr <- which(dwexac > 0)

# exact distribution using dwilcox

dw <- dwilcox(0:100, 10, 10)

# compare the two distributions:

plot(dw[tr - 10*11/2], dwexac[tr])      # should give a "perfect" line

# Ansari-Bradley

x <- rnorm(10, 2, 1)
y <- rnorm(10, 2, 2)

# exact distribution using Streitberg-Roehmel

r <- rank(c(x,y))
sc <- pmin(r, 20 - r +1)
dabexac <- dperm(0:100, sc, 10)
sum(dabexac)
tr <- which(dabexac > 0)

# exact distribution using dansari (wrapper to ansari.c in ctest)

dab <- dansari(0:100, 10, 10)

# compare the two distributions:

plot(dab[tr], dabexac[tr])

}
\keyword{exact distribution}
