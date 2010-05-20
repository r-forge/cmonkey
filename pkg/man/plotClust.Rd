\name{plotClust}
\alias{plotClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotClust(k, w.motifs = T, all.conds = F, title = NULL, o.genes = NULL, dont.plot = F, network = "all", short.names = organism == "sce", seq.type = names(mot.weights), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
~~     ~~Describe \code{k} here~~
}
  \item{w.motifs}{
~~     ~~Describe \code{w.motifs} here~~
}
  \item{all.conds}{
~~     ~~Describe \code{all.conds} here~~
}
  \item{title}{
~~     ~~Describe \code{title} here~~
}
  \item{o.genes}{
~~     ~~Describe \code{o.genes} here~~
}
  \item{dont.plot}{
~~     ~~Describe \code{dont.plot} here~~
}
  \item{network}{
~~     ~~Describe \code{network} here~~
}
  \item{short.names}{
~~     ~~Describe \code{short.names} here~~
}
  \item{seq.type}{
~~     ~~Describe \code{seq.type} here~~
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
function (k, w.motifs = T, all.conds = F, title = NULL, o.genes = NULL, 
    dont.plot = F, network = "all", short.names = organism == 
        "sce", seq.type = names(mot.weights), ...) 
{
    if (!dont.plot) 
        opar <- par(no.readonly = T)
    c <- get.clust(k)
    rows <- get.rows(k)
    if (!is.null(o.genes)) 
        rows <- c(rows, o.genes)
    if (length(rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    if (!w.motifs && !dont.plot) {
        if (all.conds) 
            plotCluster.all.conds(c, o.genes = o.genes, ...)
        else plotCluster(c, o.genes = o.genes, ...)
    }
    else {
        c$seq.type <- seq.type
        for (st in seq.type) {
            c[[st]] <- list()
            c[[st]]$motif.out <- meme.scores[[st]][[k]]
            tmp <- cluster.pclust(k, st)
            c[[st]]$e.val <- tmp$e.vals
            c[[st]]$p.clust <- tmp$p.clusts
            c[[st]]$motif.out$pssms <- lapply(c[[st]]$motif.out$meme.out, 
                "[[", "pssm")
            c[[st]]$motif.out$e.values <- c[[st]]$e.val
            if (!is.null(c[[st]]$motif.out$pv.ev)) {
                if ("gene" \%in\% colnames(c[[st]]$motif.out$pv.ev[[1]])) 
                  c[[st]]$motif.out$pv.ev[[2]] <- c[[st]]$motif.out$pv.ev[[1]]
                else if ("gene" \%in\% colnames(c[[st]]$motif.out$pv.ev[[2]])) 
                  c[[st]]$motif.out$pv.ev[[2]] <- c[[st]]$motif.out$pv.ev[[2]]
                c[[st]]$motif.out$pv.ev[[1]] <- cbind(p.value = meme.scores[[st]]$all.pv[, 
                  k], e.value = if (!is.null(meme.scores[[st]]$all.ev)) 
                  meme.scores[[st]]$all.ev[, k]
                else NA)
                c[[st]]$motif.out$p.values <- log10(c[[st]]$motif.out$pv.ev[[1]][, 
                  "p.value"])
                names(c[[st]]$motif.out$p.values) <- rownames(c[[st]]$motif.out$pv.ev[[1]])
            }
        }
    }
    if (!is.na(mot.iters[1])) {
        c$seqs <- get.sequences(rows, distance = motif.upstream.scan[[seq.type[1]]], 
            seq.type = seq.type[1], filter = F)
        tmp <- c$seqs[rows]
        if (!is.null(tmp)) 
            names(tmp) <- rows
        attr(tmp, "start.stops") <- attr(c$seqs, "start.stops")
        c$seqs <- tmp
        rm(tmp)
    }
    else c$seqs <- NULL
    if (!is.na(net.iters[1])) {
        if (network == "all") 
            network <- names(networks)
        for (i in network) {
            if (!i \%in\% names(networks)) 
                next
            tmp.net <- networks[[i]][networks[[i]]$protein1 \%in\% 
                rows & networks[[i]]$protein2 \%in\% rows, ]
            tmp.net <- cbind(tmp.net, net = rep(i, nrow(tmp.net)))
            c$network <- if (!is.null(c$network)) 
                rbind(c$network, tmp.net)
            else tmp.net
        }
    }
    c$gene.coords <- get.long.names(rows, short = short.names)
    if (!is.null(genome.info$cog.code)) 
        c$cog.code <- genome.info$cog.code[rows]
    if (!is.null(title)) 
        c$name <- title
    if (!dont.plot) {
        plotCluster.motif(c, seqs = c$seqs, p.val.shade.cutoff = 1, 
            o.genes = o.genes, ...)
        if (!"layout" \%in\% names(list(...))) 
            par(opar)
    }
    invisible(c)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
