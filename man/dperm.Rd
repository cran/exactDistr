\name{dperm}
\alias{dperm}
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

# compute the distribution of the Wilcoxon-Test using the algorithm

x <- c(0.5, 0.6, 0.7, 0.8)
y <- c(1.1, 1.2, 1.4, 1.6)
r <- rank(c(x,y))
dexact <- dperm(10:26, r, 4)	# Attention: this is the Wilcoxon NOT U
dexact == dwilcox(0:16, 4, 4) 

# exact p-value of the Wilcoxon test for a tied sample

x <- c(0.5, 0.5, 0.6, 0.6, 0.7, 0.8, 0.9)
y <- c(0.5, 1.0, 1.2, 1.2, 1.4, 1.5, 1.9, 2.0)
r <- rank(c(x,y))
pexact <- pperm(sum(r[seq(along=x)]), r, 7)
 
}
\keyword{exact distribution}
