\name{cluster.summary}
\alias{cluster.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cluster.summary(e.cutoff = 0.01, nrow.cutoff = 5, seq.type = "upstream", plot = F, sort = c("score.norm", "score", "resid", "e.value1", "e.value2", "nrow"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e.cutoff}{
~~     ~~Describe \code{e.cutoff} here~~
}
  \item{nrow.cutoff}{
~~     ~~Describe \code{nrow.cutoff} here~~
}
  \item{seq.type}{
~~     ~~Describe \code{seq.type} here~~
}
  \item{plot}{
~~     ~~Describe \code{plot} here~~
}
  \item{sort}{
~~     ~~Describe \code{sort} here~~
}
  \item{\dots}{
~~     ~~Describe \code{\dots} here~~
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
function (e.cutoff = 0.01, nrow.cutoff = 5, seq.type = "upstream", 
    plot = F, sort = c("score.norm", "score", "resid", "e.value1", 
        "e.value2", "nrow"), ...) 
{
    ms <- meme.scores[[seq.type]]
    score <- sapply(1:k.clust, function(k) mean(r.scores[get.rows(k), 
        k], na.rm = T, trim = 0.01)) * resid.scaling[iter] + 
        sapply(1:k.clust, function(k) mean(mot.scores[get.rows(k), 
            k], na.rm = T, trim = 0.01)) * mot.scaling[iter] + 
        sapply(1:k.clust, function(k) mean(net.scores[get.rows(k), 
            k], na.rm = T, trim = 0.01)) * net.scaling[iter]
    nrow <- tabulate(unlist(apply(row.membership, 1, unique)), 
        k.clust)
    out <- data.frame(k = 1:k.clust, nrow = nrow, score = score, 
        resid = sapply(1:k.clust, cluster.resid, varNorm = T), 
        consensus1 = sapply(1:k.clust, function(k) if (length(ms[[k]]) <= 
            2) 
            ""
        else pssm.to.string(ms[[k]]$meme.out[[1]]$pssm)), e.value1 = sapply(1:k.clust, 
            function(k) if (length(ms[[k]]) <= 2) 
                Inf
            else ms[[k]]$meme.out[[1]]$e.value), consensus2 = sapply(1:k.clust, 
            function(k) if (length(ms[[k]]) <= 2) 
                ""
            else if (length(ms[[k]]$meme.out) == 1) 
                ""
            else pssm.to.string(ms[[k]]$meme.out[[2]]$pssm)), 
        e.value2 = sapply(1:k.clust, function(k) if (length(ms[[k]]) <= 
            2) 
            Inf
        else if (length(ms[[k]]$meme.out) <= 1) 
            Inf
        else ms[[k]]$meme.out[[2]]$e.value))
    ss <- smooth.spline(nrow, score, spar = 0.4)
    score.norm <- score - predict(ss, nrow)$y + out$resid
    out <- cbind(out[, 1:3], score.norm, out[, 4:ncol(out)])
    if (!is.na(sort[1])) 
        out <- out[order(out[[sort[1]]]), ]
    if (!is.na(e.cutoff)) 
        out <- out[out$e.value1 <= e.cutoff | out$e.value2 <= 
            e.cutoff, ]
    if (!is.na(nrow.cutoff)) 
        out <- out[out$nrow >= nrow.cutoff, ]
    if (plot) {
        plot(out$resid, log10(-log10(out$e.value1)), typ = "n")
        text(out$resid, log10(-log10(out$e.value1)), lab = out$consensus1, 
            cex = 0.7, xpd = NA, pos = 1)
        text(out$resid, log10(-log10(out$e.value1)), lab = rownames(out), 
            cex = 0.7, xpd = NA, col = "red")
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
