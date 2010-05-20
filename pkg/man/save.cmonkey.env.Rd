\name{save.cmonkey.env}
\alias{save.cmonkey.env}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
save.cmonkey.env(env = NULL, file = NULL, restore = T, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
~~     ~~Describe \code{env} here~~
}
  \item{file}{
~~     ~~Describe \code{file} here~~
}
  \item{restore}{
~~     ~~Describe \code{restore} here~~
}
  \item{verbose}{
~~     ~~Describe \code{verbose} here~~
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
function (env = NULL, file = NULL, restore = T, verbose = T) 
{
    if (is.null(env)) {
        for (i in ls(env = globalenv())) if (is.environment(get(i, 
            globalenv())) && "cmonkey" \%in\% class(get(i, globalenv()))) 
            save.cmonkey.env(get(i, globalenv()), file, restore)
        return(invisible())
    }
    tmp <- list()
    if (is.null(file)) 
        file <- paste(env$cmonkey.filename, ".RData", sep = "")
    if (env$big.matrices > 0 && ("package:ff" \%in\% search())) {
        for (i in c("row.scores", "mot.scores", "net.scores", 
            "r.scores", "rr.scores", "col.scores", "net.scores", 
            "row.memb", "col.memb")) {
            if (!is.null(env[[i]]) && "ff" \%in\% class(env[[i]])) {
                tmp[[i]] <- env[[i]]
                env[[i]] <- env[[i]][, ]
            }
        }
        for (i in names(env$meme.scores)) {
            for (j in c("all.pv", "all.ev")) {
                if (!is.null(env$meme.scores[[i]][[j]]) && "ff" \%in\% 
                  class(env$meme.scores[[i]][[j]])) {
                  if (is.null(tmp[[i]])) 
                    tmp[[i]] <- list()
                  tmp[[i]][[j]] <- env$meme.scores[[i]][[j]]
                  env$meme.scores[[i]][[j]] <- env$meme.scores[[i]][[j]][, 
                    ]
                }
            }
        }
    }
    cat("Saving environment to", file, "\n")
    save.image(file)
    if (length(tmp) > 0 && restore && env$big.matrices > 0) {
        for (i in c("row.scores", "mot.scores", "net.scores", 
            "r.scores", "rr.scores", "col.scores", "net.scores", 
            "row.memb", "col.memb")) if (!is.null(tmp[[i]])) 
            env[[i]] <- tmp[[i]]
        for (i in names(env$meme.scores)) {
            if (is.null(tmp[[i]])) 
                next
            for (j in c("all.pv", "all.ev")) if (!is.null(tmp[[i]][[j]])) 
                env$meme.scores[[i]][[j]] <- tmp[[i]][[j]]
        }
    }
    invisible(env)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
