\name{cmonkey.init}
\alias{cmonkey.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~  ~~function to do ... ~~
}
\description{
~~  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cmonkey.init(env = NULL, ...)
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
    if (!exists("cmonkey.params")) {
        cmonkey.params <- new.env(hash = T)
    }
    if (file.exists("cmonkey-funcs.R")) {
        tmp.e <- new.env(hash = T)
        sys.source("cmonkey-funcs.R", envir = tmp.e)
    }
    else {
        tmp.e <- environment(cMonkey:::cmonkey)
    }
    if (!is.null(env) && (is.list(env) || is.environment(env))) {
        for (i in names(env)) assign(i, env[[i]])
        if (is.list(env)) 
            env <- NULL
    }
    for (i in ls(tmp.e)) {
        f2 <- NULL
        if ((!is.null(env) && exists(i, envir = env, inherit = F))) {
            f <- try(get(i, envir = env))
            f2 <- try(get(i, envir = tmp.e))
        }
        else if (exists(i, envir = .GlobalEnv, inherit = F)) {
            f <- try(get(i, envir = .GlobalEnv))
            f2 <- try(get(i, envir = tmp.e))
        }
        else if (exists(i)) {
            f <- try(get(i))
            f2 <- try(get(i, envir = tmp.e))
        }
        else {
            f <- try(get(i, envir = tmp.e))
        }
        if (class(f) == "function") 
            environment(f) <- sys.frames()[[length(sys.frames())]]
        assign(i, f)
        if (!is.null(f2) && class(f2) == "function" && object.size(f2) != 
            object.size(f)) {
            environment(f2) <- sys.frames()[[length(sys.frames())]]
            assign(paste("super", i, sep = "."), f2)
        }
    }
    rm(f, f2, tmp.e)
    if (!is.null(env)) 
        for (i in ls(env)) assign(i, get(i, env))
    args <- mget(names(formals()), env = as.environment(-1))
    for (i in names(args)) if (!i \%in\% c("...", "env")) 
        set.param(i, args[[i]])
    for (i in names(list(...))) if (i != "env") 
        set.param(i, list(...)[[i]])
    rm(args)
    if (sink.number() > 0) 
        for (i in 1:sink.number()) try(sink(), silent = T)
    set.param("save.logfile", FALSE)
    if (save.logfile != FALSE) 
        sink(save.logfile, split = T, append = (exists("dont.init") && 
            dont.init) || (exists("is.inited") && !is.inited))
    if (!exists("organism")) {
        cat("\033[31mWARNING: No organism was set; using \"hpy\".\033[0m\n")
        organism <- "hpy"
        Sys.sleep(3)
    }
    set.param("organism", organism)
    if (!exists("ratios") && exists("resid.weights")) {
        ratios <- lapply(names(resid.weights), get)
        names(ratios) <- names(resid.weights)
    }
    if ((exists("ratios") && !is.null(ratios))) {
        if (is.matrix(ratios) || is.data.frame(ratios) || is.character(ratios)) 
            ratios <- list(ratios = load.ratios(ratios))
        else ratios <- lapply(ratios, function(r) as.matrix(load.ratios(r)))
        attr(ratios, "rnames") <- sort(unique(unlist(lapply(ratios, 
            rownames))))
        attr(ratios, "cnames") <- sort(unique(unlist(lapply(ratios, 
            colnames))))
        attr(ratios, "nrow") <- length(attr(ratios, "rnames"))
        attr(ratios, "ncol") <- length(attr(ratios, "cnames"))
        for (n in names(ratios)) {
            attr(ratios[[n]], "maxRowVar") <- mean(apply(ratios[[n]], 
                1, var, use = "pair"), na.rm = T)
            attr(ratios[[n]], "all.colVars") <- apply(ratios[[n]], 
                2, var, use = "pair", na.rm = T)
        }
        rm(n)
    }
    if (!is.null(env)) 
        assign("ratios", ratios, envir = env)
    set.param("cog.org", "?")
    set.param("rsat.species", "?")
    set.param("n.iter", 3000)
    set.param("n.clust.per.row", 2)
    if (exists("ratios") && !is.null(ratios)) {
        set.param("k.clust", round(attr(ratios, "nrow") * n.clust.per.row/20))
    }
    else {
        set.param("k.clust", 100)
    }
    set.param("n.clust.per.col", if (attr(ratios, "ncol") >= 
        60) 
        round(k.clust/2)
    else round(k.clust * 2/3))
    set.param("resid.iters", seq(1, n.iter, by = 2))
    set.param("meme.iters", seq(100, n.iter, by = 100))
    set.param("mot.iters", seq(100, n.iter, by = 10))
    set.param("net.iters", seq(1, n.iter, by = 7))
    set.param("resid.scaling", 1)
    set.param("resid.weights", c(ratios = 1))
    set.param("mot.scaling", seq(0, 1, length = n.iter/2))
    set.param("mot.weights", c(upstream = 1))
    set.param("net.scaling", seq(0, 0.5, length = n.iter/2))
    set.param("net.weights", c(string = 0.5, operons = 0.5))
    set.param("grouping.weights", numeric())
    set.param("plot.iters", seq(2, n.iter, by = 25))
    set.param("post.adjust", TRUE)
    set.param("parallel.cores", TRUE)
    set.param("pareto.adjust.scalings", TRUE)
    set.param("max.changes", c(rows = 0.2, cols = 5))
    set.param("merge.cutoffs", c(n = 0.3, cor = 0.975))
    set.param("fuzzy.index", 0.75 * exp(-(1:n.iter)/(n.iter/4)))
    set.param("translation.tab", NULL)
    set.param("seed.method", c(rows = "net=string:5", cols = "best"))
    set.param("maintain.seed", NULL)
    set.param("string.version", "v8.2")
    set.param("cluster.rows.allowed", c(3, 70))
    set.param("n.motifs", c(rep(1, n.iter/3), rep(2, n.iter/3)))
    set.param("meme.cmd", "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs \%1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $compute -pal=non")
    set.param("mast.cmd", "./progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.99 -seqp -remcorr")
    set.param("dust.cmd", "./progs/dust $fname")
    set.param("operon.shift", TRUE)
    set.param("bg.order", 3)
    set.param("motif.upstream.search", c(-20, 150))
    set.param("motif.upstream.scan", c(-30, 250))
    set.param("rsat.urls", c("http://rsat.ccb.sickkids.ca/", 
        "http://rsat.ulb.ac.be/rsat/", "http://embnet.ccg.unam.mx/rsa-tools"))
    set.param("stats.iters", c(1, seq(5, n.iter, by = 5)))
    set.param("cm.script.each.iter", "cm.script.each.iter.R")
    set.param("date.run", format(Sys.time(), "\%y \%b \%d \%H:\%M:\%S"))
    set.param("cmonkey.version", cm.version)
    set.param("session.info", unlist(list(R.version, Sys.info(), 
        Sys.getenv(), sessionInfo())), quiet = T)
    set.param("time.started", date())
    if (exists("ratios") && !is.null(ratios)) {
        set.param("cmonkey.filename", paste("cmonkey", cmonkey.version, 
            organism, paste(sapply(ratios, dim), collapse = "x"), 
            gsub(" ", "_", date.run), sep = "_"))
    }
    else {
        set.param("cmonkey.filename", paste("cmonkey", cmonkey.version, 
            organism, "0x0", gsub(" ", "_", date.run), sep = "_"))
    }
    op <- options(digits.secs = 10)
    set.param("rnd.seed", as.integer(substr(gsub("[\\-\\:\\. ]", 
        "", as.character(Sys.time()), perl = T), 12, 20)))
    options(op)
    rm(op)
    set.seed(rnd.seed)
    set.param("big.matrices", 50 * 2^20)
    if (!exists("rsat.species") || rsat.species == "?" || is.na(rsat.species)) {
        err <- dlf("data/KEGG/KEGG_taxonomy.txt", "ftp://ftp.genome.jp/pub/kegg/genes/taxonomy")
        if (class(err) != "try-error") {
            tab <- read.delim("data/KEGG/KEGG_taxonomy.txt", 
                sep = "\t", comment = "#", head = F, as.is = T)
            rsat.spec <- as.character(subset(tab, V2 == organism, 
                select = "V4", drop = T))[1]
            rm(tab)
            if (any(strsplit(rsat.spec, "")[[1]] == "(")) 
                rsat.spec <- gsub("\\s\\(.*\\)", "", rsat.spec, 
                  perl = T)
        }
        rsat.spec <- gsub(" ", "_", rsat.spec, fixed = T)
        if (!file.exists("data/RSAT_genomes_listing.txt")) {
            require(RCurl)
            tmp <- strsplit(getURL(paste(rsat.urls[1], "/data/genomes/", 
                sep = "")), "\n")[[1]]
            writeLines(tmp, con = "data/RSAT_genomes_listing.txt")
        }
        vals <- character()
        if (file.exists("data/RSAT_genomes_listing.txt")) {
            tmp <- readLines("data/RSAT_genomes_listing.txt")
            vals <- grep(rsat.spec, tmp, fixed = T, val = T)
        }
        if (length(vals) <= 0) {
            cat("\033[31mCould not find correct organism for RSAT... will try to guess...\033[0m\n")
            max.dist <- 0.5
            vals <- rep("", 2)
            while (length(vals) > 1) {
                vals <- agrep(rsat.spec, tmp, ignore = T, max.dist = max.dist, 
                  val = T)
                max.dist <- max.dist - 0.01
            }
            if (length(vals) == 1) {
                rsat.spec <- strsplit(vals, "[<>/]")[[1]][8]
                cat("Found one match:", rsat.spec, "...\n")
                cat("\033[31mIf this is not correct, you're not quite out of luck -- set the 'rsat.species' parameter manually.\033[0m\n")
            }
        }
        set.param("rsat.species", rsat.spec, override = T)
        dlf(paste("data/STRING/species.", string.version, ".txt", 
            sep = ""), paste("http://string.embl.de/newstring_download/species.", 
            string.version, ".txt", sep = ""))
        rm(tmp, rsat.spec, err, vals)
    }
    else {
        set.param("rsat.species", rsat.species)
    }
    if (!exists("taxon.id") || taxon.id == "?" || is.na(taxon.id) || 
        length(taxon.id) <= 0) {
        fname <- dlf("data/GO/proteome2taxid", "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid")
        tab <- read.delim(gzfile("data/GO/proteome2taxid"), head = F)
        taxon.id <- subset(tab, V1 == gsub("_", " ", rsat.species))$V2
        if (length(taxon.id) <= 0) 
            taxon.id <- subset(tab, grepl(gsub("_", " ", rsat.species), 
                V1))$V2[1]
        set.param("taxon.id", taxon.id, override = T)
        rm(tab, fname)
    }
    if (!exists("cog.org") || cog.org == "?" || is.na(cog.org)) {
        tmp <- strsplit(organism, "")[[1]]
        tmp[1] <- toupper(tmp[1])
        cog.o <- paste(tmp, sep = "", collapse = "")
        set.param("cog.org", cog.o, override = T)
        rm(cog.o, tmp)
    }
    else {
        set.param("cog.org", cog.org)
    }
    cat("\033[32mOrganism is", organism, cog.org, rsat.species, 
        taxon.id, "\033[0m\n")
    genome.loc <- paste(rsat.urls[1], "/data/genomes/", rsat.species, 
        "/genome/", sep = "")
    fname <- paste("data/", rsat.species, "/organism.tab", sep = "")
    err <- dlf(fname, paste(genome.loc, "/organism.tab", sep = ""))
    org.tab <- readLines(fname)
    org.tab <- strsplit(org.tab[length(org.tab)], "\t")[[1]]
    is.eukaryotic <- any(grepl("Eukaryota", org.tab))
    cat("Is eukaryote:", is.eukaryotic, "\n")
    rm(err, org.tab, genome.loc, fname)
    if (is.eukaryotic) {
        cat("Organism is a eukaryote; presuming there are no operons.\n")
        set.param("is.eukaryotic", TRUE, override = T)
        set.param("operon.shift", FALSE, override = T)
        set.param("remove.low.complexity.subseqs", TRUE, override = T)
        if ("operons" \%in\% names(net.weights)) {
            net.weights <- net.weights[names(net.weights) != 
                "operons"]
            set.param("net.weights", net.weights, override = T)
        }
    }
    if (get.parallel(100, verbose = T)$mc) 
        on.exit(try(kill(children(), SIGKILL)), add = T)
    on.exit({
        if (sink.number() > 0) for (i in 1:sink.number()) try(sink(), 
            silent = T)
    })
    if (sum(net.weights, na.rm = T) > 0) 
        net.weights <- net.weights/sum(net.weights, na.rm = T)
    if (sum(resid.weights, na.rm = T) > 0) 
        resid.weights <- resid.weights/sum(resid.weights, na.rm = T)
    if (sum(mot.weights, na.rm = T) > 0) 
        mot.weights <- mot.weights/sum(mot.weights, na.rm = T)
    for (i in c("n.motifs", "meme.cmd", "mast.cmd", "operon.shift", 
        "bg.order", "motif.upstream.search", "motif.upstream.scan")) {
        v <- get(i)
        if (all(names(mot.weights) \%in\% names(v))) 
            next
        if (is.vector(v) && length(v) > 1) 
            v <- list(upstream = v)
        else names(v)[1] <- "upstream"
        for (n in names(mot.weights)[!names(mot.weights) \%in\% 
            names(v)]) {
            if (is.list(v)) 
                v[[n]] <- v[[1]]
            else if (is.vector(v)) 
                v[n] <- v[1]
            names(v) <- names(mot.weights)
        }
        assign(i, v)
    }
    rm(v)
    if (!is.null(env)) 
        for (i in ls()) {
            if (i \%in\% c("i", "env")) 
                next
            tmp <- get(i)
            if (class(tmp) == "function") 
                environment(tmp) <- env
            assign(i, tmp, envir = env)
        }
    if (!is.na(rsat.species) && (!exists("genome.info") || genome.info$species != 
        rsat.species)) {
        cat("Initializing genome info for organism", organism, 
            "\n")
        set.param("no.genome.info", FALSE)
        genome.info <- get.genome.info()
        if (!is.null(env)) 
            assign("genome.info", genome.info, envir = env)
        if (is.na(taxon.id) || length(taxon.id) <= 0) {
            taxon.id <- genome.info$taxon.id
            set.param("taxon.id", taxon.id, override = T)
        }
        genome.info$operons <- NULL
        if ((operon.shift || "operons" \%in\% names(net.weights)) && 
            !no.genome.info) {
            tmp.operons <- try(get.operon.predictions("microbes.online"))
            if (class(tmp.operons) == "try-error") {
                cat("Could not fetch operons file. Assuming it doesn't exist (eukaryote?)\n")
                set.param("is.eukaryotic", TRUE, override = T)
                set.param("operon.shift", FALSE, override = T)
                operon.shift[1:length(operon.shift)] <- FALSE
                set.param("remove.low.complexity.subseqs", TRUE, 
                  override = T)
                if ("operons" \%in\% names(net.weights)) {
                  net.weights <- net.weights[names(net.weights) != 
                    "operons"]
                  set.param("net.weights", net.weights, override = T)
                }
            }
            else {
                genome.info$operons <- tmp.operons
            }
            rm(tmp.operons)
            if (!is.null(env)) 
                assign("genome.info", genome.info, envir = env)
        }
        if (exists("genome.info") && !is.null(genome.info$feature.names) && 
            (!exists("ratios") || is.null(ratios))) 
            tmp <- toupper(subset(genome.info$feature.names, 
                type == "primary", select = "names", drop = T))
        else if (!no.genome.info) 
            tmp <- toupper(subset(genome.info$feature.names, 
                type == "primary" | names \%in\% attr(ratios, "rnames"), 
                select = "names", drop = T))
        else if (exists("ratios") && !is.null(ratios)) 
            tmp <- toupper(attr(ratios, "rnames"))
        qqq <- sapply(1:4, function(nch) max(table(substr(tmp, 
            1, nch)))/length(tmp))
        nch <- 0
        if (any(qqq > 0.6)) {
            nch <- which(qqq > 0.6)
            nch <- nch[length(nch)]
        }
        else if (any(qqq > 0.4)) {
            nch <- which(qqq > 0.4)
            nch <- nch[length(nch)]
        }
        if (nch > 0) {
            prefix <- names(which.max(table(substr(tmp, 1, nch))))
            cat("Assuming gene names have common prefix '", prefix, 
                "'.\n", sep = "")
            genome.info$gene.prefix <- prefix
        }
        else {
            cat("Could not find a common gene identifier prefix. This only matters if there's no expression matrix.\n")
            prefix <- genome.info$gene.prefix <- NA
        }
        if (!is.null(env)) 
            assign("genome.info", genome.info, envir = env)
        if (!is.na(prefix) && (!exists("ratios") || is.null(ratios))) {
            cat("WARNING: No ratios matrix -- will generate an 'empty' one with all known ORFs as probes.\n")
            rows <- unique(as.character(subset(genome.info$feature.names, 
                grepl(paste("^", prefix, sep = ""), names, ignore = T, 
                  perl = T), select = "names", drop = T)))
            ratios <- list(ratios = t(t(rep(NA, length(rows)))))
            rownames(ratios$ratios) <- rows
            attr(ratios, "rnames") <- sort(unique(rows))
            rm(rows)
            attr(ratios, "nrow") <- length(attr(ratios, "rnames"))
            attr(ratios, "ncol") <- 1
        }
        rm(nch, prefix, tmp, qqq)
        networks <- list()
        if (!is.na(net.iters) && any(net.iters \%in\% 1:n.iter)) {
            if (file.exists("data/STRING/string.csv")) {
                networks[["string"]] <- read.csv("data/STRING/string.csv", 
                  row.names = 1, header = TRUE)
            }
            else if (length(grep("string", names(net.weights))) > 
                0) {
                if ("string" \%in\% names(net.weights) || "string.combined" \%in\% 
                  names(net.weights)) {
                  if ("string.combined" \%in\% names(net.weights)) 
                    names(net.weights)[names(net.weights) == 
                      "string.combined"] <- "string"
                  string <- get.STRING.links(genome.info$org.id$V1[1], 
                    detailed = F)
                  string <- subset(string, combined_score >= 
                    500)
                  cat("Read in", nrow(string), "STRING edges that pass cutoff (500); weight =", 
                    net.weights["string"], "\n")
                  string$combined_score <- 1000 * exp(string$combined_score/1000)/exp(1)
                  networks[["string"]] <- string
                  rm(string)
                }
                if (length(grep("string.", names(net.weights))) > 
                  0) {
                  string <- get.STRING.links(genome.info$org.id$V1[1], 
                    detailed = T)
                  for (n in grep("string.", names(net.weights), 
                    val = T)) {
                    tp <- strsplit(n, ".", fixed = T)[[1]][2]
                    if (tp \%in\% colnames(string)) {
                      str <- string[, c("protein1", "protein2", 
                        tp)]
                      colnames(str)[3] <- "combined_score"
                      str <- subset(str, combined_score >= 500)
                      cat("Read in", nrow(str), n, "edges that pass cutoff (500); weight =", 
                        net.weights[n], "\n")
                      str$combined_score <- 1000 * exp(string$combined_score/1000)/exp(1)
                      networks[[n]] <- str
                    }
                  }
                  rm(string, tp, str)
                }
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if ("operons" \%in\% names(net.weights) && !is.null(genome.info$operons)) {
                cat("Converting operon predictions into a network...\n")
                tmp <- tapply(genome.info$operons$gene, genome.info$operons$head)
                names(tmp) <- genome.info$operons$gene
                mc <- get.parallel(length(unique(tmp)))
                out.sif <- do.call(rbind, mc$apply(unique(tmp), 
                  function(j) {
                    whch <- which(tmp == j)
                    gs <- names(whch)
                    if (length(gs) <= 1 || length(gs) > attr(ratios, 
                      "nrow")/20) 
                      return(NULL)
                    tmp.sif <- t(combn(gs, 2))
                    tmp.sif <- tmp.sif[tmp.sif[, 1] != tmp.sif[, 
                      2], , drop = F]
                    data.frame(protein1 = tmp.sif[, 1], protein2 = tmp.sif[, 
                      2], combined_score = rep(1000, nrow(tmp.sif)))
                  }))
                out.sif$combined_score <- rep(1000, nrow(out.sif))
                colnames(out.sif) <- c("protein1", "protein2", 
                  "combined_score")
                networks[["operons"]] <- out.sif
                rm(tmp, mc, out.sif)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (length(grep("prolinks", names(net.weights))) > 
                0) {
                prolinks.links <- get.prolinks.links(org.id = genome.info$org.id$V2[1])
                for (i in names(prolinks.links)) {
                  networks[[paste("prolinks", i, sep = ".")]] <- prolinks.links[[i]]
                  cat("Read in", nrow(prolinks.links[[i]]), i, 
                    "Prolinks edges; weight =", net.weights["prolinks"], 
                    "\n")
                }
                rm(prolinks.links, i)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (length(grep("predictome", names(net.weights))) > 
                0) {
                cat("Reading in predictome links from http://predictome.bu.edu/data/\n")
                pred.links <- get.predictome.links(org.id = organism)
                for (i in names(pred.links)) {
                  networks[[paste("pred", i, sep = ".")]] <- pred.links[[i]]
                  cat("Read in", nrow(pred.links[[i]]), i, "Predictome edges; weight =", 
                    net.weights["prolinks"], "\n")
                }
                rm(pred.links, i)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (exists("net.weights") && length(net.weights) > 
                0) {
                for (i in names(net.weights)) {
                  if (file.exists(i)) {
                    cat("Loading sif interactions from file:", 
                      i, "; weight =", net.weights[i], "\n")
                    sif <- load.sif.interactions(i)
                  }
                  else if (exists(i)) {
                    cat("Using network '", i, "' that exists in memory already; weight = ", 
                      net.weights[i], "\n", sep = "")
                    sif <- get(i)
                    if (ncol(sif) == 2) 
                      sif <- cbind(sif, rep(1, nrow(sif)))
                    colnames(sif) <- c("protein1", "protein2", 
                      "combined_score")
                  }
                  else {
                    next
                  }
                  networks[[basename(i)]] <- sif
                  rm(sif)
                }
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (exists("grouping.weights") && length(grouping.weights) > 
                0) {
                if (exists("net.weights")) 
                  net.weights <- c(net.weights, grouping.weights)
                else net.weights <- grouping.weights
                for (i in names(grouping.weights)) {
                  if (file.exists(i)) {
                    cat("Loading groupings from file:", i, "; weight =", 
                      grouping.weights[i], "\n")
                    sif <- load.sif.interactions(i)
                  }
                  else {
                    cat("Using groupings from '", i, "' that exists in memory already; weight = ", 
                      grouping.weights[i], "\n", sep = "")
                    sif <- get(i)
                    if (ncol(sif) == 2) 
                      sif <- cbind(sif, combined_score = rep(1, 
                        nrow(sif)))
                  }
                  colnames(sif) <- c("group", "protein", "combined_score")
                  if (sum(unique(as.character(sif$protein)) \%in\% 
                    attr(ratios, "rnames")) < sum(unique(as.character(sif$group)) \%in\% 
                    attr(ratios, "rnames"))) {
                    sif <- sif[, c(2, 1, 3)]
                    colnames(sif) <- c("group", "protein", "combined_score")
                  }
                  sif <- sif[order(sif$group), ]
                  tmp <- tapply(sif$protein, sif$group)
                  names(tmp) <- as.character(sif$protein)
                  cat("Converting", length(unique(tmp)), "groupings to a network (this may take a while for big grouping files)...")
                  mc <- get.parallel(length(unique(tmp)))
                  out.sif <- mc$apply(unique(tmp), function(j) {
                    whch <- which(tmp == j)
                    gs <- names(whch)
                    if (length(gs) <= 1 || length(gs) > attr(ratios, 
                      "nrow")/20) 
                      return(NULL)
                    ws <- sif$combined_score[whch]
                    names(ws) <- gs
                    tmp.sif <- t(combn(gs, 2))
                    tmp.sif <- tmp.sif[tmp.sif[, 1] != tmp.sif[, 
                      2], , drop = F]
                    tmp.sif <- data.frame(protein1 = tmp.sif[, 
                      1], protein2 = tmp.sif[, 2], combined_score = (ws[tmp.sif[, 
                      1]] + ws[tmp.sif[, 2]])/2)
                    rownames(tmp.sif) <- NULL
                    if (j\%\%100 == 0) 
                      cat(j)
                    cat(".")
                    tmp.sif
                  })
                  cat(length(unique(tmp)), "... ")
                  out.sif <- do.call(rbind, out.sif)
                  colnames(out.sif) <- c("protein1", "protein2", 
                    "combined_score")
                  networks[[basename(i)]] <- out.sif
                  cat("DONE\n")
                }
                rm(sif, tmp, out.sif, i, mc)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (exists("networks")) {
                for (n in names(networks)) {
                  nn <- networks[[n]]
                  if (nrow(nn) <= 0) {
                    cat("WARNING: no edges in network", n, "... skipping.\n")
                    if (length(grep(n, seed.method[1])) > 0) {
                      cat("\tALSO, we have to change the row seeding method from", 
                        seed.method, "to 'kmeans'.\n")
                      seed.method[1] <- "kmeans"
                      set.param("seed.method", seed.method, override = T)
                    }
                    next
                  }
                  nodes <- unique(c(as.character(nn$protein1), 
                    as.character(nn$protein2)))
                  cat(nrow(nn), "edges,", length(nodes), "nodes in network", 
                    n, "\n")
                  nn <- subset(nn, as.character(protein1) != 
                    as.character(protein2))
                  dupes <- duplicated(nn[, c("protein1", "protein2")])
                  if (sum(dupes) > 0) {
                    cat("Merging", sum(dupes), "duplicate edges in network", 
                      n, "; this could take a while for networks with lots of nodes...\n")
                    tmp.nn <- subset(nn, dupes)
                    dupe.nodes <- unique(c(as.character(tmp.nn$protein1), 
                      as.character(tmp.nn$protein2)))
                    if (length(dupe.nodes) < 6000) {
                      tmp <- tapply(tmp.nn$combined_score, tmp.nn[, 
                        c("protein1", "protein2")], sum, na.rm = T)
                      tmp2 <- which(!is.na(tmp), arr = T)
                      nn.new <- data.frame(protein1 = rownames(tmp)[tmp2[, 
                        1]], protein2 = colnames(tmp)[tmp2[, 
                        2]], combined_score = tmp[tmp2])
                      rm(tmp, tmp2)
                      nn <- rbind(nn.new, nn)
                      rm(nn.new)
                      nn <- nn[!duplicated(nn[, c("protein1", 
                        "protein2")]), ]
                    }
                    rm(tmp.nn, dupe.nodes)
                  }
                  if (exists("ratios") && !is.null(ratios) && 
                    !any(nodes \%in\% attr(ratios, "rnames"))) {
                    if (median(nchar(nodes)) > median(nchar(attr(ratios, 
                      "rnames"))) && any(substr(nodes, 1, median(nchar(attr(ratios, 
                      "rnames")))) \%in\% attr(ratios, "rnames"))) {
                      nn$protein1 <- substr(as.character(nn$protein1), 
                        1, median(nchar(attr(ratios, "rnames"))))
                      nn$protein2 <- substr(as.character(nn$protein2), 
                        1, median(nchar(attr(ratios, "rnames"))))
                      nodes <- unique(c(as.character(nn$protein1), 
                        as.character(nn$protein2)))
                    }
                    if (!is.null(genome.info$synonyms)) {
                      rr <- attr(ratios, "rnames")[!attr(ratios, 
                        "rnames") \%in\% nodes]
                      if (length(rr) > 0) {
                        cat("Reconciling network", n, length(rr), 
                          "node names with probe names...\n")
                        syns <- get.synonyms(rr)
                        mc <- get.parallel(length(syns))
                        is.there <- unlist(mc$apply(syns, function(i) any(i \%in\% 
                          nodes)))
                        syns <- syns[is.there]
                        nnc1 <- as.character(nn$protein1)
                        nnc2 <- as.character(nn$protein2)
                        nnc1.t <- !nnc1 \%in\% attr(ratios, "rnames")
                        nnc2.t <- !nnc2 \%in\% attr(ratios, "rnames")
                        mc <- get.parallel(2)
                        tmp <- mc$apply(1:2, function(ii) {
                          for (i in names(syns)) {
                            if (ii == 1) 
                              nnc1[nnc1.t & nnc1 \%in\% syns[[i]]] <- i
                            else nnc2[nnc2.t & nnc2 \%in\% syns[[i]]] <- i
                          }
                          if (ii == 1) 
                            return(nnc1)
                          else return(nnc2)
                        })
                        nnc1 <- tmp[[1]]
                        nnc2 <- tmp[[2]]
                        rm(tmp, nnc1.t, nnc2.t)
                        cat(sum(!is.there), "probes have no nodes in", 
                          n, "network (but", sum(attr(ratios, 
                            "rnames") \%in\% nodes, na.rm = T) + 
                            sum(is.there), "do)\n")
                        nn$protein1 <- nnc1
                        nn$protein2 <- nnc2
                        tmp <- nnc1 \%in\% attr(ratios, "rnames") & 
                          nnc2 \%in\% attr(ratios, "rnames")
                        nn <- subset(nn, tmp == TRUE)
                        rm(tmp, syns, is.there, nnc1, nnc2, nnc1.t, 
                          nnc2.t, tmp, rr, i)
                      }
                    }
                  }
                  else {
                    cat(sum(!attr(ratios, "rnames") \%in\% nodes), 
                      "probes have no nodes in", n, "network (but", 
                      sum(attr(ratios, "rnames") \%in\% nodes, 
                        na.rm = T), "do)\n")
                  }
                  ttmp <- nn[, c(2, 1, 3)]
                  colnames(ttmp) <- colnames(nn)
                  nn <- rbind(nn, ttmp)
                  rm(ttmp)
                  nn <- nn[!duplicated(nn[, c("protein1", "protein2")]), 
                    ]
                  cat(n, "network filtered, symmetrized and uniquified:", 
                    nrow(nn), "edges.\n")
                  networks[[n]] <- nn
                  if (!is.null(env)) 
                    assign("networks", networks, envir = env)
                }
                rm(n, nn, nodes, dupes)
                if (length(networks) > 1) {
                  sums <- sapply(networks, function(n) sum(n$combined_score, 
                    na.rm = T))
                  ms <- min(sums[sums > 0], na.rm = T)
                  if (length(sums) > 0 && !is.na(ms)) 
                    for (n in names(networks)) networks[[n]]$combined_score <- networks[[n]]$combined_score/sums[n] * 
                      ms
                  rm(n, sums, ms)
                }
                if (!is.null(env)) 
                  assign("networks", networks, envir = env)
            }
            names(net.weights) <- basename(names(net.weights))
        }
        if (!is.null(genome.info$genome.seqs)) {
            genome.info$all.upstream.seqs <- genome.info$bg.list <- list()
            for (i in names(mot.weights)) {
                cat("Pre-computing all", i, "seqs for background distribution (", 
                  motif.upstream.scan[[i]], ")...\n")
                genome.info$all.upstream.seqs[[i]] <- get.sequences(attr(ratios, 
                  "rnames"), seq.type = i, distance = motif.upstream.scan[[i]])
                if (!is.null(env)) 
                  assign("genome.info", genome.info, envir = env)
                cat(sum(!attr(ratios, "rnames") \%in\% names(genome.info$all.upstream.seqs[[i]])), 
                  "probes have no", i, "sequence.\n")
                if (!is.na(bg.order[i])) {
                  cat("Pre-computing", i, "residue bg distrib (order=", 
                    bg.order[i], ")...\n")
                  tmp.seqs <- if (!is.null(genome.info$all.upstream.seqs[[i]])) 
                    genome.info$all.upstream.seqs[[i]]
                  else get.sequences(attr(ratios, "rnames"), 
                    distance = motif.upstream.search[[i]], seq.type = i)
                  capture.output(genome.info$bg.list[[i]] <- mkBgFile(tmp.seqs, 
                    order = bg.order[i], use.rev.comp = grepl("-revcomp", 
                      meme.cmd[i])))
                  rm(tmp.seqs)
                }
                else {
                  cat("NOT USING a global sequence background distribution!\n")
                }
                if (!is.null(env)) 
                  assign("genome.info", genome.info, envir = env)
            }
        }
        if (!no.genome.info) {
            cat("Loading COG functional codes (for plotting), org. code", 
                cog.org, ": trying NCBI whog file...\n")
            genome.info$cog.code <- get.COG.code(cog.org)
        }
        cat(sum(!is.na(genome.info$cog.code)), "genes have a COG code (", 
            sum(is.na(genome.info$cog.code)), "do not)\n")
    }
    iter <- 0
    meme.scores <- clusterStack <- list()
    for (i in names(mot.weights)) {
        meme.scores[[i]] <- list()
        meme.scores[[i]][[k.clust + 1]] <- ""
    }
    stats <- row.scores <- col.scores <- mot.scores <- net.scores <- r.scores <- NULL
    old.row.membership <- old.col.membership <- r.scores <- NULL
    if (!exists("favorite.cluster")) 
        favorite.cluster <- function() min(which(tabulate(row.membership) > 
            cluster.rows.allowed[1] * 2))
    resid.scaling <- extend.vec(resid.scaling)
    mot.scaling <- extend.vec(mot.scaling)
    net.scaling <- extend.vec(net.scaling)
    n.motifs <- lapply(n.motifs, extend.vec)
    fuzzy.index <- extend.vec(fuzzy.index)
    is.inited <- TRUE
    if (is.null(env)) 
        env <- new.env(hash = T, parent = globalenv())
    else parent.env(env) <- globalenv()
    parent.env(cmonkey.params) <- env
    attr(env, "class") <- c("environment", "cmonkey")
    for (i in ls()) {
        if (i \%in\% c("i", "env")) 
            next
        tmp <- get(i)
        if (class(tmp) == "function") 
            environment(tmp) <- env
        assign(i, tmp, envir = env)
    }
    if (exists("favorite.cluster")) 
        env$favorite.cluster <- favorite.cluster
    environment(env$favorite.cluster) <- env
    if (exists("cm.func.each.iter")) {
        env$cm.func.each.iter <- cm.func.each.iter
        environment(env$cm.func.each.iter) <- env
        try(env$cm.func.each.iter())
    }
    cat("INITIALIZATION IS COMPLETE.\n")
    env$iter <- env$iter + 1
    if (env$n.iter == 3000) 
        env$n.iter <- 2000
    invisible(env)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
