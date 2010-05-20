\name{plot.stats}
\alias{plot.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.stats(iter = stats$iter[nrow(stats)], plot.clust = NA, new.dev = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iter}{
~~     ~~Describe \code{iter} here~~
}
  \item{plot.clust}{
~~     ~~Describe \code{plot.clust} here~~
}
  \item{new.dev}{
~~     ~~Describe \code{new.dev} here~~
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
function (iter = stats$iter[nrow(stats)], plot.clust = NA, new.dev = T, 
    ...) 
{
    row.memb <- t(apply(row.membership, 1, function(i) 1:k.clust \%in\% 
        i))
    opar <- par(no.readonly = T)
    tmp.scale <- round(1/mean(row.memb, na.rm = T)/4)
    if (new.dev) {
        if (length(dev.list()) < 1) 
            dev.new()
        dev.set(2)
    }
    layout(matrix(c(1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8, 
        10, 7, 9, 11), byrow = T, ncol = 3))
    par(mar = c(3, 3, 2, 0.1), mgp = c(3, 1, 0) * 0.5)
    stats <- stats[stats[, "iter"] <= iter, , drop = F]
    try(matplot(stats[, "iter"], stats[, grep("resid", colnames(stats), 
        val = T)], typ = "l", xlab = "iter", ylab = "Mean resid", 
        main = sprintf("Iter: \%d", iter)), silent = T)
    if ((nn <- length(grep("resid", colnames(stats)))) > 1) 
        legend("bottomleft", legend = gsub("resid.", "", grep("resid", 
            colnames(stats), val = T)), lwd = 1, bty = "n", col = 1:nn, 
            lty = 1:nn, cex = 0.5)
    rs <- row.scores[]
    rs[rs < -220] <- min(rs[rs > -220], na.rm = T)
    h <- try(hist(rs, breaks = 50, main = NULL, xlab = "Ratios scores"))
    if (class(h) != "try-error") {
        try(hist(rep(rs[row.memb], tmp.scale), breaks = h$breaks, 
            col = "red", border = "red", add = T), silent = T)
        try(hist(rs, breaks = h$breaks, add = T), silent = T)
    }
    if (!is.null(mot.scores)) {
        ms <- mot.scores[]
        ms[ms < -20] <- min(ms[ms > -20], na.rm = T)
        ms[ms >= 0] <- NA
        h <- try(hist(ms, breaks = 50, main = NULL, xlab = "Motif scores"))
        if (class(h) != "try-error") {
            try(hist(rep(ms[row.memb], tmp.scale * 3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ms, breaks = h$breaks, add = T), silent = T)
        }
        try(matplot(stats[, "iter"], stats[, grep("p.clust", 
            colnames(stats), val = T)], typ = "l", xlab = "iter", 
            ylab = "Mean motif p-value", main = sprintf("Motif scaling: \%.3f", 
                mot.scaling[max(1, iter - 1)])), silent = T)
        if ((nn <- length(grep("p.clust", colnames(stats)))) > 
            1) 
            legend("bottomleft", legend = gsub("p.clust.", "", 
                grep("p.clust", colnames(stats), val = T)), lwd = 1, 
                bty = "n", col = 1:nn, lty = 1:nn, cex = 0.5)
    }
    if (!is.null(net.scores)) {
        ns <- net.scores[]
        ns[ns < -20] <- min(ns[ns > -20], na.rm = T)
        ns[ns >= 0] <- NA
        ns[, ] <- -log10(-ns)
        tmp.scale <- ceiling(tmp.scale * mean(!is.na(ns), na.rm = T))
        h <- try(hist(ns, breaks = 50, main = NULL, xlab = "-log10(-Network scores)"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ns[row.memb], tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ns, breaks = h$breaks, add = T))
        }
        try(matplot(stats[, "iter"], stats[, grep("net.", colnames(stats), 
            val = T, fixed = T)], typ = "l", xlab = "iter", ylab = "Mean net-score", 
            main = sprintf("Net scaling: \%.3f", net.scaling[max(1, 
                iter - 1)])), silent = T)
        if ((nn <- length(grep("net.", colnames(stats)))) > 1) 
            try(legend("bottomleft", legend = gsub("net.", "", 
                grep("net.", colnames(stats), val = T)), lwd = 1, 
                bty = "n", col = 1:nn, lty = 1:nn, cex = 0.5), 
                silent = T)
    }
    clusterStack <- get.clusterStack(ks = 1:k.clust)
    resids <- sapply(clusterStack, "[[", "resid")
    try(hist(resids[resids <= 1.5], main = NULL, xlab = "Cluster Residuals", 
        xlim = c(0, 1.5), breaks = k.clust/4), silent = T)
    if (!is.null(mot.scores)) {
        plot.all.clusterMotifPositions <- function(ks = 1:k.clust, 
            mots = 1, e.cutoff = 1, p.cutoff = 0.05, seq.type = "upstream", 
            breaks = 100, ...) {
            ms <- meme.scores[[seq.type]]
            posns <- as.vector(unlist(sapply(ms[ks], function(i) i$pv.ev[[1]]$posns)))
            pvals <- as.vector(unlist(sapply(ms[ks], function(i) i$pv.ev[[1]]$pvals)))
            imots <- as.vector(unlist(sapply(ms[ks], function(i) i$pv.ev[[1]]$mots)))
            clusts <- as.vector(unlist(sapply(ms[ks], function(i) rep(i$k, 
                if (is.null(i$pv.ev[[1]])) 0 else nrow(i$pv.ev[[1]])))))
            evals <- sapply(1:length(imots), function(i) ms[[clusts[i]]]$meme.out[[abs(imots[i])]]$e.value)
            psns <- posns[evals < e.cutoff & pvals < p.cutoff & 
                abs(imots) \%in\% mots] - motif.upstream.scan[[seq.type]][2]
            h <- hist(psns, breaks = breaks, ...)
            dd <- density(psns, bw = 5)
            lines(dd$x, dd$y * max(h$counts)/max(dd$y) * 0.9, 
                col = "red")
            invisible(data.frame(clusts, posns, pvals, imots, 
                evals))
        }
        plot.all.clusterMotifPositions(xlab = "Position upstream", 
            main = "Positions of motif #1", ...)
        if (any(sapply(e$meme.scores$upstream[1:e$k.clust], function(i) length(i$meme.out)) == 
            2)) 
            plot.all.clusterMotifPositions(mots = 2, xlab = "Position upstream", 
                main = "Positions of motif #2", ...)
    }
    n.rows <- tabulate(unlist(apply(row.membership, 1, unique)))
    try(hist(n.rows, main = NULL, xlab = "Cluster Nrows", breaks = k.clust/4, 
        xlim = c(-5, max(n.rows, na.rm = T))), silent = T)
    n.cols <- tabulate(unlist(apply(col.membership, 1, unique)))
    try(hist(n.cols, main = NULL, xlab = "Cluster Ncols", breaks = k.clust/4, 
        xlim = c(-5, attr(ratios, "ncol"))), silent = T)
    if (!is.na(plot.clust)) {
        if (new.dev) {
            if (length(dev.list()) < 2) 
                dev.new()
            dev.set(3)
        }
        try(plotClust(plot.clust, T, cex = 0.7), silent = T)
        if (new.dev) {
            if (length(dev.list()) < 3) 
                dev.new()
            dev.set(4)
        }
        try(plot.scores(plot.clust), silent = T)
    }
    par(opar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
