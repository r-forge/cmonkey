\name{clusters.w.func}
\alias{clusters.w.func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clusters.w.func(func, ks = 1:k.clust, short = F, max.rows = 999, p.val = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
~~     ~~Describe \code{func} here~~
}
  \item{ks}{
~~     ~~Describe \code{ks} here~~
}
  \item{short}{
~~     ~~Describe \code{short} here~~
}
  \item{max.rows}{
~~     ~~Describe \code{max.rows} here~~
}
  \item{p.val}{
~~     ~~Describe \code{p.val} here~~
}
}
\details{
~~  ~~ If necessary, more details than the description above ~~
}
\value{
~~  ~Describe the value returned
~~  If it is a LIST, use
~~  \item{comp1 }{Description of 'comp1'}
~~  \item{comp2 }{Description of 'comp2'}
~~ ...
}
\references{
~~ ~put references to the literature/web site here ~
}
\author{
~~  ~~who you are~~
}
\note{
~~  ~~further notes~~
}

~~ ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
~~ ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (func, ks = 1:k.clust, short = F, max.rows = 999, p.val = F) 
{
    if (p.val) 
        n2 <- length(grep(func, long.names, perl = T, ignore.case = T))
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        if (!p.val) {
            if (length(get.rows(i)) >= max.rows) 
                NA
            else length(grep(func, get.long.names(get.rows(i), 
                short = short), perl = T, ignore.case = T))
        }
        else {
            phyper(length(grep(func, get.long.names(get.rows(i), 
                short = short), perl = T, ignore.case = T)), 
                n2, attr(ratios, "nrow") - n2, length(get.rows(i)), 
                lower = F)/length(ks)
        }
    }))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
