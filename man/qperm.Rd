\name{qperm}
\alias{qperm}
\title{Distribution of Permutation Tests}
\description{
Density, distribution function and quantile function 
for the distribution of permutation tests using the Shift-Algorithm by
Streitberg and R\"ohmel.
}
\usage{
dperm(x, scores, m, paired=NULL, tol = 0.01, fact=NULL)
pperm(q, scores, m, paired=NULL, tol = 0.01, fact=NULL)
qperm(p, scores, m, paired=NULL, tol = 0.01, fact=NULL)
rperm(n, scores, m)
}
\arguments{
  \item{x, q}{ vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{scores}{ranks, midranks or real valued scores of the observations 
of the \code{x} (first \code{m} elements) and \code{y} sample. }
  \item{m}{ sample size of the \code{x} sample. If \code{m = length(x)}
scores of paired observations are assumed.}
  \item{paired}{ logical. Indicates if paired observations are used. Only
needed for the distribution of the total sum of the scores (which has mass 1
at the point \code{sum(scores)})}.
  \item{tol}{ real. real valued scores are mapped into integers by multiplication.
Make sure that the absolute difference between the "true" quantile and the
approximated quantile is less than \code{tol}. This might not be possible
due to memory/time limitations.} 
  \item{fact}{ real. If \code{fact} is given, real valued scores are mapped into
integers using \code{fact} as factor. \code{tol} is ignored.}
  \item{n}{number of observations.}
}
\details{
  The exact distribution of the sum of the first \code{m} scores is
evaluated using the Shift-Algorithm by Streitberg and R\"ohmel under the
hypothesis of exchangeability.

The algorithm is defined for positive integer valued scores only. 
There are two ways dealing with real valued scores. 
First, one can try to find integer
valued scores that lead to quantiles which differ not more than \code{tol}
from the quantiles induced by the original scores. This can be done as
follows.  

Without loss of generality let \eqn{a_i > 0} denote real valued scores and
\eqn{f} a positive factor. Let \eqn{R_i = a_i - round(f \cdot a_i)}.  
Then 

\deqn{ \sum_{i=1}^m f \cdot a_i = \sum_{i=1}^m round(f \cdot a_i) - R_i. }

Clearly, the maximum difference between \eqn{\sum_{i=1}^m f \cdot a_i} and
\eqn{\sum_{i=1}^n round(f \cdot a_i)} is given by \eqn{|\sum_{i=1}^m R_{(i)}|} or
\eqn{|\sum_{i=m+1}^N R_{(i)}|}, respectively. Therfore one searches for
\eqn{f} with 

\deqn{ \max(|\sum_{i=1}^m R_{(i)}|, |\sum_{i=m+1}^N R_{(i)}|)/f \le tol.}

If \eqn{f} induces more that 20.000 columns in the Streitberg-R\"ohmel
Shift-Algorithm, \eqn{f} is restricted to the largest integer that does not. 

The second idea is to map the scores itself into
\eqn{\{1, \dots, N\}}. This induces additional ties, but the shape of the
scores is very similar. That means we do not try to approximate something
but use a different test (with integer scores), serving for the same purpose
(the shape of the scores is very similar). 

}
\value{
\code{dperm} gives the density, \code{pperm} gives the distribution
function and \code{qperm} gives the quantile function. \code{rperm} is
just a one-line wrapper to \code{\link{sample}}.

}
\references{ 

	Bernd Streitberg and Joachim R\"ohmel (1986). 
	Exact Distributions For Permutations and Rank Tests: 
	An Introduction to Some Recently Published Algorithms. 
	\emph{Statistical Software Newsletter} \bold{12}, No. 1, 10--17.

	Bernd Streitberg and Joachim R\"ohmel (1987)
	Exakte Verteilungen f\"ur Rang- und Randomisierungstests 
	im allgemeinen $c$-Stichprobenfall.
	\emph{EDV in Medizin und Biologie} \bold{18}, No. 1, 12--19.
}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{
\dontrun{

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
su
if (su != 1) stop("sum(dwexac) not equal 1")

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

# real scores are allowed (but only result in an approximation)
# e.g. v.d. Waerden test

x <- rnorm(10)
y <- rnorm(10)
N <- length(x) + length(y)
r <- rank(c(x,y))
scores <- qnorm(r/(N+1))
X <- sum(scores[seq(along=x)])  # <- v.d. Waerden normal quantile statistic

# critical value, two-sided test

abs(qperm(0.025, scores, length(x)))

# p-values

p1 <- pperm(abs(X), scores, length(x))
min(p1, 1-p1)

# generate integer valued scores with the same shape as normal quantile
# scores, this no longer v.d.Waerden, but something very similar

N <- 20
scores <- scores - min(scores)
scores <- round(scores*N/max(scores))

X <- sum(scores[seq(along=x)])
p2 <- pperm(abs(X), scores, length(x))
min(p2, 1-p2)

# compare p1 and p2

p1 - p2

# check if [qp]wilcox and [pq]ansari are equal with [qp]perm

source(system.file("test1.R", pkg="exactDistr")[1])
source(system.file("test2.R", pkg="exactDistr")[1])

}
}
\keyword{exact distribution}
