\name{update.cmonkey.env}
\alias{update.cmonkey.env}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update.cmonkey.env(env)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
~~     ~~Describe \code{env} here~~
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
function (env) 
{
    if (file.exists("cmonkey-funcs.R")) {
        tmp.e <- new.env()
        sys.source("cmonkey-funcs.R", envir = tmp.e)
        sys.source("cmonkey-postproc.R", envir = tmp.e)
    }
    else {
        tmp.e <- environment(cMonkey:::cmonkey)
    }
    for (i in ls(tmp.e)) {
        if (i \%in\% c("DATE", "VERSION")) 
            next
        f <- try(get(i, envir = tmp.e))
        f2 <- try(get(paste("super", i, sep = "."), envir = env), 
            silent = T)
        if (class(f) == "function") {
            environment(f) <- env
            if (class(f2) != "function") 
                assign(i, f)
            else assign(paste("super", i, sep = "."), f)
        }
    }
    rm(f, f2, tmp.e, i)
    for (i in ls()) {
        if (i \%in\% c("i", "env")) 
            next
        f <- get(i)
        if (is.function(f)) 
            assign(i, f, env)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
