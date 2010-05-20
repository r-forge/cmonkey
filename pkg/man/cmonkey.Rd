\name{cmonkey}
\alias{cmonkey}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cmonkey(env = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{env}{
~~     ~~Describe \code{env} here~~
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
function (env = NULL, ...) 
{
    if (((is.null(list(...)$dont.init) || !list(...)$dont.init) && 
        (is.null(env$dont.init) || !env$dont.init) && (!exists("dont.init") || 
        !dont.init)) || is.null(env) || is.null(env$genome.info)) {
        env <- cmonkey.init(env, ...)
    }
    else {
        if (sink.number() > 0) 
            for (i in 1:sink.number()) try(sink(), silent = T)
        if (env$save.logfile != FALSE) 
            sink(env$save.logfile, split = T, append = T)
    }
    cat("\033[31mTIME STARTED:", env$time.started, "\033[0m\n")
    if ((!exists("row.membership", envir = env) || all(env$row.membership == 
        0) || nrow(env$row.membership) != attr(env$ratios, "nrow") || 
        nrow(env$col.membership) != attr(env$ratios, "ncol")) && 
        exists("ratios", envir = env)) {
        if (!exists("rnd.seed", envir = env$cmonkey.params)) {
            op <- options(digits.secs = 10)
            tmp.rnd.seed <- as.integer(substr(gsub("[\\-\\:\\. ]", 
                "", as.character(Sys.time()), perl = T), 12, 
                20))
            options(op)
            rm(op)
            cat("RESETTING RANDOM SEED: ")
            env$set.param("rnd.seed", tmp.rnd.seed, env$cmonkey.params)
            env$rnd.seed <- env$cmonkey.params$rnd.seed
            set.seed(env$rnd.seed)
            rm(tmp.rnd.seed)
        }
        if (!is.null(env$ratios) && attr(env$ratios, "ncol") > 
            1) {
            cat("Seeding all clusters using methods:", env$seed.method, 
                "\n")
            tmp <- env$seed.clusters(env$k.clust, seed.method = env$seed.method["rows"], 
                col.method = env$seed.method["cols"])
        }
        else {
            cat("Seeding all clusters using methods: rnd rnd\n")
            tmp <- env$seed.clusters(env$k.clust, seed.method = "rnd", 
                col.method = "rnd")
        }
        env$row.membership <- tmp$row.membership
        env$col.membership <- tmp$col.membership
        rm(tmp)
    }
    iter <- env$iter
    while (iter <= env$n.iter) {
        env$iter <- iter
        env$cmonkey.one.iter(env)
        iter <- iter + 1
    }
    if (!is.na(env$plot.iters) && (iter \%in\% env$plot.iters || 
        (iter - 1) \%in\% env$plot.iters)) 
        try(env$plot.stats(iter, plot.clust = env$favorite.cluster()))
    env$iter <- iter <- iter - 1
    print(env$cluster.summary())
    if (env$post.adjust == TRUE) {
        env$pre.adjusted.row.membership <- env$row.membership
        env$adjust.all.clusters(env)
        if (env$big.matrices > 0) 
            env$ffify.env(env)
        gc()
    }
    parent.env(env) <- globalenv()
    parent.env(env$cmonkey.params) <- env
    print(env$cluster.summary())
    env$set.param("time.ended", date(), env$cmonkey.params)
    env$time.ended <- date()
    cat("\033[31mTIME ENDED:", env$time.ended, "\033[0m\n")
    invisible(env)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
