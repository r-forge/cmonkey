DATE <-
"Mon Aug 22 11:36:42 2011"
VERSION <-
"5.0.0b"
.onLoad <-
function (libname, pkgname) 
{
    packageStartupMessage("Loading ", pkgname, " version ", VERSION, 
        " (", DATE, ")")
    packageStartupMessage("Copyright (C) David J Reiss, Institute for Systems Biology; dreiss@systemsbiology.org.")
    packageStartupMessage("http://baliga.systemsbiology.net/cmonkey")
    vers <- try(readLines("http://baliga.systemsbiology.net/cmonkey/VERSION"), 
        silent = T)
    if (class(vers) != "try-error") {
        vers <- gsub(" ", "", vers)
        if (vers != VERSION) 
            packageStartupMessage("\nYou are not using the most current version of cMonkey.\nPlease consider upgrading to v", 
                vers, " via:\n\n> download.file( \"http://baliga.systemsbiology.net/cmonkey/cMonkey_", 
                vers, ".tar.gz\", \n\t\t\"cMonkey_", vers, ".tar.gz\" )\n> install.packages( \"cMonkey_", 
                vers, ".tar.gz\", repos=NULL )\n\nOr by following the instructions on the cMonkey website.")
        else packageStartupMessage("Congratulations! You are using the latest version of cMonkey.\n")
    }
    else {
        packageStartupMessage("Could not check to see if you are using the latest version of cMonkey.")
    }
}
DEBUG <-
function (...) 
{
    message(...)
}
adjust.all.clusters <-
function (env, ks = 1:env$k.clust, force.motif = T, ...) 
{
    old.stats <- env$stats
    tmp <- env$row.col.membership.from.clusterStack(env$clusterStack)
    row.membership <- tmp$r
    col.membership <- tmp$c
    mc <- env$get.parallel(length(ks))
    new.rm <- mc$apply(ks, function(k) env$adjust.clust(k, row.membership, 
        ...)$r)
    rm <- do.call(cbind, new.rm)
    for (i in 1:nrow(rm)) {
        tmp <- unique(rm[i, rm[i, ] != 0])
        rm[i, ] <- c(tmp, rep(0, ncol(rm) - length(tmp)))
    }
    rm <- rm[, apply(rm, 2, sum) != 0, drop = F]
    colnames(rm) <- NULL
    env$clusterStack <- lapply(1:env$k.clust, function(k) list(rows = rownames(which(rm == 
        k, arr = T)), cols = env$clusterStack[[k]]$cols))
    env$clusterStack <- env$get.clusterStack(ks = 1:k.clust)
    env$post.adjust <- FALSE
    env$cmonkey.one.iter(env, dont.update = T, force.row = T, 
        force.col = T, force.motif = if (force.motif & !no.genome.info) 
            "run.meme", force.net = T)
    print(rbind(OLD = old.stats[nrow(old.stats), ], NEW = env$stats[nrow(env$stats), 
        ]))
    invisible(env)
}
adjust.clust <-
function (k, row.memb = get("row.membership"), expand.only = T, 
    limit = 100, scores = "r.scores", quant.cutoff = 0.33, force.expand = 0) 
{
    if (scores == "rr.scores" || scores == "r.scores") {
        tmp <- get.combined.scores(quant = T)
        r.scores <- tmp$r
        if (scores == "rr.scores") {
            scores <- get.density.scores(ks = 1:k.clust)$r
            scores <- 1 - scores[, ]
        }
        else {
            scores <- r.scores
        }
        rm(r.scores)
    }
    else {
        scores <- get(scores)
    }
    get.rows2 <- function(k, rm) rownames(which(rm == k, arr = T))
    scores <- scores[, ]
    old.rows <- get.rows(k)
    if (force.expand == 0) {
        wh <- names(which(scores[which(!attr(ratios, "rnames") %in% 
            old.rows), k] < quantile(scores[old.rows, k], quant.cutoff, 
            na.rm = T)))
    }
    else {
        expand.only <- TRUE
        wh <- names(sort(scores[!attr(ratios, "rnames") %in% 
            old.rows, k], decreasing = F)[1:force.expand])
    }
    if (length(wh) > limit) {
        warning("Surpassing limit.")
        return(invisible(list(r = row.memb)))
    }
    else if (length(wh) <= 0) 
        return(invisible(list(r = row.memb)))
    tries <- 0
    while (length(wh) > 0 && tries < 50) {
        wh2 <- names(which.max(scores[wh, k]))
        wh2.scores <- scores[wh2, row.memb[wh2, ]]
        wh2a <- names(which.max(scores[get.rows2(k, rm = row.memb), 
            k]))
        for (col in 1:ncol(row.memb)) if (all(row.memb[wh2, col] == 
            0)) 
            break
        if (col == ncol(row.memb) && any(row.memb[wh2, col] != 
            0)) {
            row.memb <- cbind(row.memb, rep(0, nrow(row.memb)))
            col <- col + 1
        }
        row.memb[wh2, col] <- k
        if (!expand.only) 
            row.memb[wh2a, row.memb[wh2a, ] == k] <- 0
        if (force.expand == 0) {
            wh <- names(which(scores[which(!attr(ratios, "rnames") %in% 
                get.rows2(k, rm = row.memb)), k] < quantile(scores[get.rows2(k, 
                rm = row.memb), k], quant.cutoff, na.rm = T)))
        }
        else {
            wh <- wh[!wh %in% wh2]
        }
        if (length(get.rows2(k, rm = row.memb)) > cluster.rows.allowed[2]) 
            break
        tries <- tries + 1
    }
    new.rows <- get.rows2(k, rm = row.memb)
    if (any(!new.rows %in% old.rows) || any(!old.rows %in% new.rows)) 
        cat("ADJUSTED CLUSTER:", k, length(old.rows), length(new.rows), 
            sum(!old.rows %in% new.rows), "\n")
    row.memb <- t(apply(row.memb, 1, function(i) c(i[i != 0], 
        i[i == 0])))
    row.memb <- row.memb[, apply(row.memb, 2, sum) != 0, drop = F]
    colnames(row.memb) <- NULL
    invisible(list(r = row.memb))
}
agglom <-
function (src = "MOTC_1", srcType = "motif.cluster", targetType = "gene", 
    path = "motif,bicluster", p.val = T, q.val = "BH", cond.filter = NULL, 
    cond.filter.frac = 0.5, verbose = F) 
{
    by <- unlist(strsplit(path, ","))
    orig.src <- src
    orig.srcType <- srcType
    for (b in by) {
        text <- paste("bys <- unlist( get.", b, "s( ", srcType, 
            "=src ) )", sep = "")
        if (verbose) 
            print(text)
        eval(parse(text = text))
        if (!is.null(cond.filter) && b == "bicluster") {
            if (verbose) 
                print("cc <- get.conditions( biclust=bys )")
            cc <- get.conditions(biclust = bys)
            tmp <- sapply(cc, function(i) mean(cond.filter %in% 
                i) >= cond.filter.frac)
            if (!any(tmp)) 
                warning("No biclusters pass the cond.filter! Perhaps too many conditions?")
            else bys <- bys[tmp]
        }
        src <- bys
        srcType <- b
    }
    text <- paste("out <- get.", targetType, "s( ", by[length(by)], 
        "=bys )", sep = "")
    if (verbose) 
        print(text)
    eval(parse(text = text))
    if (targetType == "tf") 
        out <- lapply(out, lapply, names)
    tab <- sort(table(unlist(out)), decreasing = T)
    out2 <- t(t(tab))
    colnames(out2) <- "count"
    tab2 <- tab[tab > 1]
    pvals <- qvals <- numeric()
    if (p.val && length(tab2) > 0 && by[length(by)] == "bicluster") {
        bics1 <- as.list(names(tab))
        names(bics1) <- bics1
        if (targetType != "bicluster") {
            text <- paste("bics1 <- get.biclusters( ", targetType, 
                "=names( tab2 ) )", sep = "")
            if (verbose) 
                print(text)
            eval(parse(text = text))
        }
        n.bics <- e$k.clust
        if (!is.null(cond.filter)) {
            bics <- table(unlist(get.biclusters(cond = cond.filter)))
            bics <- bics[bics >= length(cond.filter) * cond.filter.frac]
            bics1 <- lapply(bics1, function(i) i[i %in% names(bics)])
            n.bics <- length(bics)
        }
        pvals <- phyper(tab2, length(out), n.bics, sapply(bics1, 
            length)[names(tab2)], lower = F)
        if (!is.na(q.val)) {
            if (q.val == TRUE) 
                q.val <- "BH"
            pv <- c(pvals, rep(1, nrow(out2) - length(pvals)))
            qvals <- p.adjust(pv, q.val)
        }
    }
    if (p.val) {
        out2 <- cbind(out2, p.value = c(pvals, rep(1, nrow(out2) - 
            length(pvals))))
        if (!is.na(q.val)) 
            out2 <- cbind(out2, q.value = c(qvals, rep(1, nrow(out2) - 
                length(qvals))))
    }
    out2 <- as.data.frame(out2)
    attr(out2, "total") <- length(bys)
    out2
}
agglom2 <-
function (src = "MOTC_1", srcType = get.type(src), srcPath = "motif,bicluster", 
    target = "MOTC_170", targetType = get.type(target), targetPath = "motif,bicluster", 
    cond.filter = NULL, cond.filter.frac = 0.5, verbose = F) 
{
    by <- unlist(strsplit(srcPath, ","))
    orig.src <- src
    orig.srcType <- srcType
    for (b in by) {
        text <- paste("bys1 <- unlist( get.", b, "s( ", srcType, 
            "=src ) )", sep = "")
        if (verbose) 
            print(text)
        eval(parse(text = text))
        if (!is.null(cond.filter) && b == "bicluster") {
            if (verbose) 
                print("cc <- get.conditions( biclust=bys1 )")
            cc <- get.conditions(biclust = bys1)
            tmp <- sapply(cc, function(i) mean(cond.filter %in% 
                i) >= cond.filter.frac)
            if (!any(tmp)) 
                warning("No biclusters pass the cond.filter! Perhaps too many conditions?")
            else bys1 <- bys1[tmp]
        }
        src <- bys1
        srcType <- b
    }
    by <- unlist(strsplit(targetPath, ","))
    orig.target <- target
    orig.targetType <- targetType
    for (b in by) {
        text <- paste("bys2 <- unlist( get.", b, "s( ", targetType, 
            "=target ) )", sep = "")
        if (verbose) 
            print(text)
        eval(parse(text = text))
        if (!is.null(cond.filter) && b == "bicluster") {
            if (verbose) 
                print("cc <- get.conditions( biclust=bys2 )")
            cc <- get.conditions(biclust = bys2)
            tmp <- sapply(cc, function(i) mean(cond.filter %in% 
                i) >= cond.filter.frac)
            if (!any(tmp)) 
                warning("No biclusters pass the cond.filter! Perhaps too many conditions?")
            else bys2 <- bys2[tmp]
        }
        target <- bys2
        targetType <- b
    }
    out <- data.frame(src = orig.src, target = orig.target, nSrc = length(bys1), 
        nTarget = length(bys2), nBoth = sum(bys1 %in% bys2), 
        p.value = phyper(sum(bys1 %in% bys2), length(bys1), e$k.clust - 
            length(bys1), length(bys2), lower = F))
    out
}
all.dna.seqs <-
function (l, lett = c("G", "A", "T", "C"), as.matrix = F) 
{
    n.lett <- length(lett)
    out <- sapply(1:l, function(ll) rep(as.vector(sapply(lett, 
        function(i) rep(i, n.lett^(ll - 1)))), n.lett^(l - ll)))
    if (as.matrix) 
        return(out)
    apply(out, 1, paste, collapse = "")
}
blast.align <-
function (seqs1, seqs2, addl.params = "", full.out = F, verbose = F, 
    unlink = T, bl2seq.cmd = sprintf("%s/bl2seq", progs.dir), 
    formatdb.cmd = sprintf("%s/formatdb", progs.dir), blast.cmd = sprintf("%s/blastall", 
        progs.dir)) 
{
    file1 <- my.tempfile("blast.")
    file2 <- my.tempfile("blast.")
    tmp.log.file <- my.tempfile("formatdb.log")
    if (length(seqs1) == 1 && length(seqs2) == 1) {
        cat(">seq1\n", seqs1, "\n", sep = "", file = file1)
        cat(">seq2\n", seqs2, "\n", sep = "", file = file2)
        cmd <- sprintf("%s -i %s -j %s -p blastn %s", bl2seq.cmd, 
            file1, file2, addl.params)
        if (!full.out) 
            cmd <- paste(cmd, "-D 1")
        if (verbose) 
            cat(cmd, "\n")
        output <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    }
    else {
        if (length(seqs1) == 1) 
            cat(">seq1\n", seqs1, "\n", sep = "", file = file1)
        else cat(paste(">", names(seqs1), "\n", seqs1, sep = ""), 
            file = file1, sep = "\n")
        if (length(seqs2) == 1) 
            cat(">seq2\n", seqs2, "\n", sep = "", file = file2)
        else cat(paste(">", names(seqs2), "\n", seqs2, sep = ""), 
            file = file2, sep = "\n")
        if (length(seqs1) >= length(seqs2)) {
            cmd <- sprintf("%s -l %s -i %s -p F -o T", formatdb.cmd, 
                tmp.log.file, file1)
            if (verbose) 
                cat(cmd, "\n")
            system(cmd)
            cmd <- sprintf("%s -p blastn -d %s -i %s %s", blast.cmd, 
                file1, file2, addl.params)
        }
        else if (length(seqs2) > length(seqs1)) {
            cmd <- sprintf("%s -l %s -i %s -p F -o T", formatdb.cmd, 
                tmp.log.file, file2)
            if (verbose) 
                cat(cmd, "\n")
            system(cmd)
            cmd <- sprintf("%s -p blastn -d %s -i %s %s", blast.cmd, 
                file2, file1, addl.params)
        }
        if (!full.out) 
            cmd <- paste(cmd, "-m 8")
        if (verbose) 
            cat(cmd, "\n")
        output <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    }
    if (unlink) 
        try(unlink(c(paste(file1, "*", sep = ""), paste(file2, 
            "*", sep = ""), tmp.log.file)))
    return(output)
}
blast.match.seqs <-
function (seqs, match = NULL, e.cutoff = 1) 
{
    if (is.null(match)) 
        match <- seqs
    out <- parse.blast.out(blast.align(seqs, match, paste("-e", 
        e.cutoff)))
    out <- subset(out, as.character(Query.id) != as.character(Subject.id))
    out <- subset(out, as.character(Query.id) %in% names(seqs))
    out <- subset(out, as.character(Subject.id) %in% names(match))
    return(out[order(out$e.value), ])
}
cluster.loglik <-
function (k) 
{
    l.in <- rr.scores[get.rows(k), k]
    l.out <- rr.scores[!rownames(rr.scores) %in% get.rows(k), 
        k]
    sum(log(c(l.in, 1 - l.out) + 1e-99), na.rm = T)
}
cluster.meme.motif.lines <-
function (k, seq.type = names(meme.scores)[1], logodds = F) 
{
    meme.let <- c("A", "C", "G", "T")
    lines <- c("ALPHABET= ACGT", "", "strands: + -", "", "Background letter frequencies (from dataset with add-one prior applied):")
    lines <- c(lines, paste(names(unlist(genome.info$bg.list[[seq.type]][meme.let])), 
        sprintf("%.3f", unlist(genome.info$bg.list[[seq.type]][meme.let])), 
        collapse = " "))
    memeOut <- meme.scores[[seq.type]][[k]]$meme.out
    if (is.null(memeOut)) 
        return(lines)
    for (i in 1:length(memeOut)) {
        pssm <- memeOut[[i]]$pssm
        mat.type <- "letter-probability matrix"
        if (logodds) 
            mat.type <- "log-odds matrix"
        lines <- c(lines, "", sprintf("MOTIF bic_%03d_%02d", 
            k, i), sprintf("BL   MOTIF bic_%03d_%02d width=0 seqs=0", 
            k, i), sprintf("%s: alength= 4 w= %d nsites= %d E= %.3e", 
            mat.type, nrow(pssm), memeOut[[i]]$sites, memeOut[[i]]$e.value))
        if (!logodds) 
            lines <- c(lines, apply(pssm, 1, function(i) sprintf("%5.3f %5.3f %5.3f %5.3f", 
                i[1], i[2], i[3], i[4])))
        else lines <- c(lines, apply(round(log(pssm + 0.01)), 
            1, function(i) sprintf("%6d %6d %6d %6d", i[1], i[2], 
                i[3], i[4])))
    }
    lines
}
cluster.pclust <-
function (k, mot.inds = "COMBINED") 
{
    inds <- mot.inds
    if (is.null(inds)) 
        return(list(p.clusts = NA, e.vals = NA))
    if (mot.inds[1] == "COMBINED") 
        inds <- names(get("mot.weights"))
    if (is.null(inds)) 
        return(list(p.clusts = NA, e.vals = NA))
    rows <- get.rows(k)
    p.clusts <- sapply(inds, function(n) {
        ms <- meme.scores[[n]][[k]]
        if (length(rows) > 0 && !is.null(ms$pv.ev) && !is.null(ms$pv.ev[[1]])) 
            mean(log10(ms$pv.ev[[1]][rownames(ms$pv.ev[[1]]) %in% 
                rows, "p.value"]), na.rm = T)
        else NA
    })
    e.vals <- sapply(inds, function(n) {
        ms <- meme.scores[[n]][[k]]
        sapply(1:length(ms$meme.out), function(i) if (length(rows) > 
            0 && !is.null(ms$meme.out) && !is.null(ms$meme.out[[i]])) 
            ms$meme.out[[i]]$e.value
        else NA)
    })
    if (!is.matrix(e.vals)) 
        e.vals <- t(t(e.vals))
    if (mot.inds[1] == "COMBINED") {
        p.clusts <- weighted.mean(p.clusts, mot.weights[inds], 
            na.rm = T)
        e.vals <- apply(e.vals, 2, weighted.mean, mot.weights[inds], 
            na.rm = T)
    }
    else if (mot.inds[1] != "COMBINED") {
        if (length(p.clusts) < length(inds) && all(is.na(p.clusts))) 
            p.clusts <- rep(NA, length(inds))
        e.vals <- apply(e.vals, 1, function(i) if (length(i) < 
            length(inds) && all(is.na(i))) 
            rep(NA, length(inds))
        else i)
    }
    if (!is.matrix(e.vals)) 
        e.vals <- t(t(e.vals))
    if (is.matrix(e.vals) && ncol(e.vals) != length(inds)) 
        e.vals <- t(e.vals)
    if (mot.inds[1] != "COMBINED") 
        names(p.clusts) <- colnames(e.vals) <- inds
    else e.vals <- as.vector(e.vals)
    list(p.clusts = p.clusts, e.vals = e.vals)
}
cluster.resid <-
function (k, rats.inds = "COMBINED", varNorm = F, in.cols = T, 
    ...) 
{
    residual.submatrix <- function(rats, rows, cols, varNorm = F, 
        ...) {
        rows <- rows[rows %in% rownames(rats)]
        cols <- cols[cols %in% colnames(rats)]
        if (length(rows) <= 1 || length(cols) <= 1) 
            return(1)
        maxRowVar <- attr(rats, "maxRowVar")
        rats <- rats[rows, cols]
        if (is.vector(rats) || any(dim(rats) <= 1) || mean(is.na(rats)) > 
            0.95) 
            return(1)
        d.rows <- rowMeans(rats, na.rm = T)
        d.cols <- colMeans(rats, na.rm = T)
        d.all <- mean(d.rows, na.rm = T)
        rats[, ] <- rats[, ] + d.all - outer(d.rows, d.cols, 
            "+")
        average.r <- mean(abs(rats), na.rm = TRUE)
        if (varNorm && !is.null(maxRowVar)) {
            row.var <- mean(apply(rats, 1, var, use = "pairwise.complete.obs"), 
                na.rm = T)
            if (is.na(row.var) || row.var > maxRowVar) 
                row.var <- maxRowVar
            average.r <- average.r/row.var
        }
        average.r
    }
    inds <- rats.inds
    if (rats.inds[1] == "COMBINED") 
        inds <- names(get("row.weights"))
    rows <- get.rows(k)
    cols <- get.cols(k)
    resids <- sapply(ratios[inds], function(rn) {
        if (row.score.func == "orig") {
            if (in.cols) 
                residual.submatrix(rn, rows, cols, varNorm = varNorm)
            else residual.submatrix(rn, get.rows(k), colnames(rn)[!colnames(rn) %in% 
                cols], varNorm = varNorm)
        }
        else {
            if (in.cols) 
                mean(get.row.scores(k, for.rows = rows, ratios = rn, 
                  method = row.score.func))
            else mean(get.row.scores(k, cols = cols, for.rows = rows, 
                ratios = rn, method = row.score.func))
        }
    })
    if (rats.inds[1] == "COMBINED") 
        resids <- weighted.mean(resids, row.weights[inds], na.rm = T)
    if (rats.inds[1] != "COMBINED" && length(resids) < length(inds) && 
        all(is.na(resids))) {
        resids <- rep(NA, length(inds))
        names(resids) <- inds
    }
    resids
}
cluster.summary <-
function (e.cutoff = 0.01, nrow.cutoff = 5, seq.type = names(mot.weights)[1], 
    plot = F, sort = c("score.norm", "score", "resid", "e.value1", 
        "e.value2", "nrow")) 
{
    ms <- NULL
    if (!is.null(seq.type)) 
        ms <- meme.scores[[seq.type]]
    if (is.null(ms)) 
        e.cutoff <- NA
    score <- sapply(1:k.clust, function(k) mean(row.scores[get.rows(k), 
        k], na.rm = T, trim = 0.01)) * row.scaling[iter] + if (!is.null(mot.scores)) 
        sapply(1:k.clust, function(k) mean(mot.scores[get.rows(k), 
            k], na.rm = T, trim = 0.01)) * mot.scaling[iter]
    else 0 + if (!is.null(net.scores)) 
        sapply(1:k.clust, function(k) mean(net.scores[get.rows(k), 
            k], na.rm = T, trim = 0.01)) * net.scaling[iter]
    else 0
    nrow <- sapply(1:k.clust, function(k) length(get.rows(k)))
    out <- data.frame(k = 1:k.clust, nrow = nrow, score = score, 
        resid = sapply(1:k.clust, cluster.resid, varNorm = F), 
        consensus1 = sapply(1:k.clust, function(k) if (is.null(ms) || 
            is.null(ms[[k]]$meme.out) || length(ms[[k]]) <= 3) 
            ""
        else pssm.to.string(ms[[k]]$meme.out[[1]]$pssm)), e.value1 = sapply(1:k.clust, 
            function(k) if (is.null(ms) || is.null(ms[[k]]$meme.out) || 
                length(ms[[k]]) <= 3) 
                Inf
            else ms[[k]]$meme.out[[1]]$e.value), consensus2 = sapply(1:k.clust, 
            function(k) if (is.null(ms) || is.null(ms[[k]]$meme.out) || 
                length(ms[[k]]) <= 3) 
                ""
            else if (length(ms[[k]]$meme.out) == 1) 
                ""
            else pssm.to.string(ms[[k]]$meme.out[[2]]$pssm)), 
        e.value2 = sapply(1:k.clust, function(k) if (is.null(ms) || 
            is.null(ms[[k]]$meme.out) || length(ms[[k]]) <= 3) 
            Inf
        else if (length(ms[[k]]$meme.out) <= 1) 
            Inf
        else ms[[k]]$meme.out[[2]]$e.value))
    if (all(out$consensus2 == "")) 
        out$consensus2 <- out$e.value2 <- NULL
    if (!is.na(sort[1]) && sort[1] %in% colnames(out)) 
        out <- out[order(out[[sort[1]]]), ]
    if (!all(is.na(score))) {
        ss <- smooth.spline(nrow[!is.na(score)], score[!is.na(score)], 
            spar = 0.6)
        score.norm <- score - predict(ss, nrow)$y + out$resid
        out <- cbind(out[, 1:3], score.norm, out[, 4:ncol(out)])
        if (!is.na(e.cutoff)) {
            if ("e.value2" %in% colnames(out)) 
                out <- out[out$e.value1 <= e.cutoff | out$e.value2 <= 
                  e.cutoff, ]
            else out <- out[out$e.value1 <= e.cutoff, ]
        }
        if (!is.na(nrow.cutoff)) 
            out <- out[out$nrow >= nrow.cutoff, ]
        if (plot) {
            plot(out$resid, log10(-log10(out$e.value1)), typ = "n")
            text(out$resid, log10(-log10(out$e.value1)), lab = out$consensus1, 
                cex = 0.7, xpd = NA, pos = 1)
            text(out$resid, log10(-log10(out$e.value1)), lab = rownames(out), 
                cex = 0.7, xpd = NA, col = "red")
        }
    }
    out
}
cluster.tomtom.results <-
function (tt.out, seq.type = names(mot.weights)[1], e.value.cutoff = Inf, 
    p.value.cutoff = 0.05, resid.cutoff = Inf, n.cutoff = 3, 
    make.pssms = T, min.size = 3, k.cut = 0.5, return.aligned.pssms = F, 
    n.gene.weight = F, ...) 
{
    if (is.na(e.value.cutoff)) 
        e.value.cutoff <- Inf
    if (is.na(resid.cutoff)) 
        resid.cutoff <- Inf
    if (is.na(p.value.cutoff)) 
        p.value.cutoff <- Inf
    if (!is.infinite(e.value.cutoff) || !is.infinite(resid.cutoff) || 
        !is.infinite(p.value.cutoff)) 
        tt.out <- subset(tt.out, e.value1 <= e.value.cutoff & 
            e.value2 <= e.value.cutoff & resid1 <= resid.cutoff & 
            resid2 <= resid.cutoff & p.value <= p.value.cutoff)
    mot.names <- c(paste(tt.out$biclust1, tt.out$motif1, sep = "_"), 
        paste(tt.out$biclust2, abs(tt.out$motif2), sep = "_"))
    mot.tab <- sort(table(mot.names), decreasing = T)
    mot.names <- names(mot.tab)
    mot.names.2 <- cbind(paste(tt.out$biclust1, tt.out$motif1, 
        sep = "_"), paste(tt.out$biclust2, abs(tt.out$motif2), 
        sep = "_"))
    tmp <- matrix(1, nrow = length(mot.names), ncol = length(mot.names))
    rownames(tmp) <- colnames(tmp) <- mot.names
    tmp[mot.names.2] <- tt.out$p.value
    tmp[cbind(mot.names.2[, 2], mot.names.2[, 1])] <- tt.out$p.value
    tmp[tmp > p.value.cutoff] <- 1
    hc <- hclust(as.dist(tmp), method = "complete")
    if (k.cut < 1) 
        cs <- cutree(hc, h = k.cut)
    else cs <- cutree(hc, k = k.cut)
    rm(tmp)
    gc()
    cat("HERE:", range(as.integer(names(sort(table(cs), decreasing = T)))), 
        "\n")
    meme.let <- c("A", "C", "G", "T")
    mc <- get.parallel(length(table(cs)))
    out.tt <- mc$apply(as.integer(names(sort(table(cs), decreasing = T))), 
        function(i) {
            wh <- t(apply(do.call(rbind, strsplit(names(which(cs == 
                i)), "_")), 1, as.integer))
            if (nrow(wh) < n.cutoff) 
                return(NULL)
            cat(i, nrow(wh), "\n")
            wh.tmp <- apply(wh, 1, paste, collapse = "_")
            b1 <- paste(tt.out$biclust1, abs(tt.out$motif1), 
                sep = "_")
            b2 <- paste(tt.out$biclust2, abs(tt.out$motif2), 
                sep = "_")
            tt <- subset(tt.out, b1 %in% wh.tmp & b2 %in% wh.tmp)
            if (nrow(tt) <= 0) 
                return(NULL)
            if ("q.value" %in% names(tt)) 
                tt <- tt[order(tt$p.value, tt$q.value), ]
            else tt <- tt[order(tt$p.value), ]
            b1 <- paste(tt$biclust1, abs(tt$motif1), sep = "_")
            b2 <- paste(tt$biclust2, abs(tt$motif2), sep = "_")
            if (make.pssms) {
                tto <- tt
                mot.names <- unique(c(b1, b2))
                tmp1 <- as.integer(strsplit(mot.names[1], "_")[[1]])
                tmp.mot <- meme.scores[[seq.type]][[tmp1[1]]]$meme.out[[tmp1[2]]]
                pssm <- orig.pssm <- tmp.mot$pssm
                colnames(pssm) <- meme.let
                if (n.gene.weight) 
                  pssm <- pssm * tmp.mot$sites
                if (return.aligned.pssms) 
                  aligned.pssms <- list()
                if (length(mot.names) >= min.size) {
                  if ("width1" %in% names(tto)) 
                    max.width <- max(c(tto$width1, tto$width2))
                  else max.width <- max(nchar(c(as.character(tto$consensus1), 
                    as.character(tto$consensus2))))
                  for (jj in 1:max.width) pssm <- rbind(rep(0, 
                    4), pssm, rep(0, 4))
                  if (return.aligned.pssms) 
                    aligned.pssms[[mot.names[1]]] <- pssm
                  first.ind <- which(apply(pssm, 1, function(j) any(j != 
                    0)))[1]
                  orig.width <- nrow(orig.pssm)
                  tttt <- subset(tto, (biclust1 == tmp1[1] & 
                    motif1 == tmp1[2]) | (biclust2 == tmp1[1] & 
                    abs(motif2) == tmp1[2]))
                  rm(tto)
                  for (m in mot.names[2:length(mot.names)]) {
                    tmp <- as.integer(strsplit(m, "_")[[1]])
                    ttt <- unique(subset(tttt, (biclust1 == tmp[1] & 
                      abs(motif1) == tmp[2]) | (biclust2 == tmp[1] & 
                      abs(motif2) == tmp[2])))
                    if (nrow(ttt) <= 0) 
                      next
                    else if (nrow(ttt) > 1) 
                      ttt <- ttt[1, ]
                    tmp.mot <- meme.scores[[seq.type]][[tmp[1]]]$meme.out[[tmp[2]]]
                    pssm2 <- tmp.mot$pssm
                    if (n.gene.weight) 
                      pssm2 <- pssm2 * tmp.mot$sites
                    if (("orientation" %in% colnames(ttt) && 
                      ttt$orientation == "-") || ttt$motif2 < 
                      0) {
                      if (ttt$biclust2 == tmp1[1] && abs(ttt$motif2) == 
                        tmp1[2]) 
                        offset <- first.ind + orig.width - ttt$offset - 
                          nrow(pssm2)
                      else offset <- first.ind - ttt$offset
                      pssm2 <- pssm2[, 4:1][nrow(pssm2):1, ]
                    }
                    else {
                      if (ttt$biclust1 == tmp1[1] && abs(ttt$motif1) == 
                        tmp1[2]) 
                        offset <- first.ind - ttt$offset
                      else offset <- first.ind + ttt$offset
                    }
                    pssm[offset:(offset + nrow(pssm2) - 1), ] <- pssm[offset:(offset + 
                      nrow(pssm2) - 1), ] + pssm2
                    if (return.aligned.pssms) {
                      tmp.pssm <- pssm * 0
                      tmp.pssm[offset:(offset + nrow(pssm2) - 
                        1), ] <- tmp.pssm[offset:(offset + nrow(pssm2) - 
                        1), ] + pssm2
                      aligned.pssms[[m]] <- tmp.pssm
                    }
                  }
                }
                first.ind2 <- which(apply(pssm, 1, function(j) any(j != 
                  0)))[1]
                pssm <- pssm[-(1:(first.ind2 - 1)), ]
                last.ind2 <- which(apply(pssm, 1, function(j) all(j == 
                  0)))[1]
                if (!is.na(last.ind2)) 
                  pssm <- pssm[-((last.ind2 - 1):nrow(pssm)), 
                    ]
                pssm <- (pssm + 1e-09)/(max(pssm, na.rm = T) + 
                  1e-09)
                attr(tt, "combined.pssm") <- pssm
                if (return.aligned.pssms) {
                  for (m in names(aligned.pssms)) {
                    pssm <- aligned.pssms[[m]]
                    pssm <- pssm[-(1:(first.ind2 - 1)), ]
                    pssm <- pssm[-((last.ind2 - 1):nrow(pssm)), 
                      ]
                    aligned.pssms[[m]] <- pssm
                  }
                  attr(tt, "aligned.pssms") <- aligned.pssms
                }
                attr(tt, "wh") <- wh
                attr(tt, "mot.names") <- mot.names
            }
            tt
        })
    attr(out.tt, "hc.out") <- hc
    out.tt
}
clusters.w.conds <-
function (conds, ks = 1:k.clust, p.val = F) 
{
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        cols <- get.cols(i)
        if (!p.val) 
            sum(cols %in% conds)
        else phyper(sum(cols %in% conds), length(conds), attr(ratios, 
            "ncol") - length(conds), length(cols), lower = F) * 
            length(ks)
    }))
}
clusters.w.func <-
function (func, ks = 1:k.clust, short = F, max.rows = 999, p.val = F) 
{
    if (p.val) {
        long.names <- get.long.names(attr(ratios, "rnames"), 
            short = short)
        n2 <- length(grep(func, long.names, perl = T, ignore.case = T))
    }
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        rows <- get.rows(i)
        if (length(rows) <= 1) 
            return(NA)
        rows.l <- get.long.names(rows, short = short)
        if (!p.val) {
            if (length(rows) >= max.rows) 
                NA
            else length(grep(func, rows.l, perl = T, ignore.case = T))
        }
        else {
            phyper(length(grep(func, rows.l, perl = T, ignore.case = T)), 
                n2, attr(ratios, "nrow") - n2, length(rows), 
                lower = F) * length(ks)
        }
    }))
}
clusters.w.genes <-
function (genes, ks = 1:k.clust, p.val = F) 
{
    mc <- get.parallel(length(ks))
    unlist(mc$apply(ks, function(i) {
        rows <- get.rows(i)
        if (length(rows) <= 1) 
            return(NA)
        if (!p.val) 
            sum(rows %in% genes)
        else phyper(sum(rows %in% genes), length(genes), attr(ratios, 
            "nrow") - length(genes), length(rows), lower = F) * 
            length(ks)
    }))
}
cm.version <-
"5.0.0"
cmonkey <-
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
    if ((!exists("clusterStack", envir = env) || length(env$clusterStack) < 
        env$k.clust) && exists("ratios", envir = env)) 
        env$cmonkey.re.seed(env)
    while (env$iter <= env$n.iter) {
        iter <- env$iter
        env$cmonkey.one.iter(env)
    }
    if (!is.na(env$plot.iters) && (iter %in% env$plot.iters || 
        (iter - 1) %in% env$plot.iters)) 
        try(env$plotStats(iter, plot.clust = env$favorite.cluster()))
    env$iter <- iter <- env$iter - 1
    print(env$cluster.summary())
    parent.env(env) <- globalenv()
    parent.env(env$cmonkey.params) <- env
    env$clusterStack <- env$get.clusterStack(ks = 1:env$k.clust, 
        force = T)
    print(env$cluster.summary())
    env$set.param("time.ended", date(), env$cmonkey.params)
    env$time.ended <- date()
    cat("\033[31mTIME ENDED:", env$time.ended, "\033[0m\n")
    invisible(env)
}
cmonkey.init <-
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
        for (i in names(env)) if (!i %in% names(list(...))) 
            assign(i, env[[i]])
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
            if (!is.null(env)) {
                assign(paste("super", i, sep = "."), f2, envir = env)
                environment(env[[paste("super", i, sep = ".")]]) <- env
            }
        }
    }
    rm(f, f2, tmp.e)
    if (!is.null(env)) 
        for (i in ls(env)) assign(i, get(i, env))
    args <- mget(names(formals()), env = as.environment(-1))
    for (i in names(args)) if (!i %in% c("...", "env")) 
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
        message("WARNING: No organism was set; using \"hpy\".")
        organism <- "hpy"
        Sys.sleep(3)
    }
    set.param("organism", organism)
    if (!exists("ratios") && exists("row.weights")) {
        try({
            ratios <- lapply(names(row.weights), get)
            names(ratios) <- names(row.weights)
        })
    }
    if ((exists("ratios") && !is.null(ratios))) {
        if (is.matrix(ratios) || is.data.frame(ratios)) 
            ratios <- list(ratios = load.ratios(ratios))
        else if (is.character(ratios)) 
            ratios <- lapply(ratios, load.ratios)
        else ratios <- lapply(ratios, function(r) as.matrix(load.ratios(r)))
        ratios <- ratios[sapply(ratios, function(r) all(dim(r) > 
            0))]
        attr(ratios, "rnames") <- sort(unique(unlist(lapply(ratios, 
            rownames))))
        attr(ratios, "cnames") <- sort(unique(unlist(lapply(ratios, 
            colnames))))
        attr(ratios, "nrow") <- length(attr(ratios, "rnames"))
        attr(ratios, "ncol") <- length(attr(ratios, "cnames"))
        if (is.null(names(ratios))) {
            names(ratios) <- paste("ratios", 1:length(ratios), 
                sep = ".")
            if (!all(names(ratios) %in% names(row.weights))) 
                for (i in names(ratios)) row.weights[i] <- row.weights[1]
        }
        for (n in names(ratios)) {
            if (ncol(ratios[[n]]) > 1) {
                attr(ratios[[n]], "maxRowVar") <- mean(apply(ratios[[n]][, 
                  ], 1, var, use = "pair"), na.rm = T)
                attr(ratios[[n]], "all.colVars") <- apply(ratios[[n]][, 
                  ], 2, var, use = "pair", na.rm = T)
            }
        }
        rm(n)
    }
    if (exists("ratios") && is.null(names(ratios))) 
        names(ratios) <- paste("ratios", 1:length(ratios), sep = ".")
    if (!is.null(env) && exists("ratios")) 
        assign("ratios", ratios, envir = env)
    set.param("cog.org", "?")
    set.param("rsat.species", "?")
    set.param("n.iter", 2000)
    set.param("n.clust.per.row", 2)
    if (exists("ratios") && !is.null(ratios)) {
        set.param("k.clust", round(attr(ratios, "nrow") * n.clust.per.row/20))
    }
    else {
        set.param("k.clust", 100)
    }
    set.param("n.clust.per.col", if (exists("ratios") && attr(ratios, 
        "ncol") >= 60) 
        round(k.clust/2)
    else round(k.clust * 2/3))
    set.param("row.iters", seq(1, n.iter, by = 2))
    set.param("col.iters", seq(1, n.iter, by = 5))
    set.param("meme.iters", c(seq(600, 1200, by = 100), seq(1250, 
        1500, by = 50), seq(1525, 1800, by = 25), seq(1810, max(n.iter, 
        1820), by = 10)))
    set.param("mot.iters", seq(601, max(n.iter, 605), by = 3))
    set.param("net.iters", seq(1, n.iter, by = 7))
    set.param("row.scaling", 6)
    set.param("row.weights", c(ratios = 1))
    set.param("row.score.func", "orig")
    set.param("col.score.func", "orig")
    set.param("mot.scaling", seq(0, 1, length = n.iter * 3/4))
    set.param("mot.weights", c(`upstream meme` = 1))
    set.param("net.scaling", seq(0, 0.5, length = n.iter * 3/4))
    set.param("net.weights", c(string = 0.5, operons = 0.5))
    set.param("grouping.weights", numeric())
    set.param("plot.iters", seq(2, n.iter, by = 25))
    set.param("post.adjust", TRUE)
    set.param("parallel.cores", TRUE)
    set.param("parallel.cores.motif", TRUE)
    set.param("max.changes", c(rows = 0.5, cols = 5))
    set.param("cluster.rows.allowed", c(3, 70))
    set.param("merge.cutoffs", c(n = 0.3, cor = 0.975))
    set.param("fuzzy.index", 0.7 * exp(-(1:n.iter)/(n.iter/3)) + 
        0.05)
    set.param("translation.tab", NULL)
    set.param("seed.method", c(rows = "kmeans", cols = "best"))
    set.param("maintain.seed", NULL)
    set.param("n.motifs", c(rep(1, n.iter/2), rep(2, n.iter/4), 
        3))
    if (file.exists("./progs")) {
        set.param("progs.dir", "./progs/")
    }
    else if ("package:cMonkey" %in% search() && file.exists(sprintf("%s/progs/", 
        system.file(package = "cMonkey")))) {
        set.param("progs.dir", sprintf("%s/progs/", system.file(package = "cMonkey")))
    }
    else if (any(mot.scaling > 0) && (!exists("no.genome.info") || 
        !no.genome.info)) {
        message("WARNING: You do not have meme/mast/dust installed.\nTrying to install it now.\n")
        install.binaries()
        set.param("progs.dir", sprintf("%s/progs/", system.file(package = "cMonkey")))
        if ("package:cMonkey" %in% search() && !file.exists(sprintf("%s/progs/", 
            system.file(package = "cMonkey")))) 
            message("WARNING: Could not install meme. Please see the website for installation instructions.")
    }
    set.param("meme.cmd", paste(progs.dir, "meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $none", 
        sep = "/"))
    set.param("mast.cmd", sprintf("%s/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.99 -seqp -remcorr", 
        progs.dir))
    set.param("weeder.cmd", "./weederlauncher.out %s %s %s S T%d")
    set.param("spacer.cmd", c("java -Xmx1000M -jar SPACER.jar -b %s -o %s %s", 
        "java -Xmx1000M -jar SPACER.jar -l %s -o %s %s"))
    set.param("dust.cmd", sprintf("%s/dust $fname", progs.dir))
    set.param("operon.shift", TRUE)
    set.param("bg.order", 3)
    set.param("recalc.bg", TRUE)
    set.param("motif.upstream.search", c(-20, 150))
    set.param("motif.upstream.scan", c(-30, 250))
    set.param("discard.genome", TRUE)
    set.param("pareto.adjust.scalings", FALSE)
    set.param("rsat.urls", c("http://rsat.ccb.sickkids.ca/", 
        "http://rsat.ulb.ac.be/rsat/", "http://embnet.ccg.unam.mx/rsa-tools"))
    set.param("stats.iters", c(1, seq(5, n.iter, by = 5)))
    set.param("cm.script.each.iter", "cm.script.each.iter.R")
    set.param("date.run", format(Sys.time(), "%y %b %d %H:%M:%S"))
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
    set.param("rnd.seed", as.integer(substr(gsub("[-:. ]", "", 
        as.character(Sys.time())), 12, 20)))
    options(op)
    rm(op)
    set.seed(rnd.seed)
    set.param("big.memory", FALSE)
    set.param("big.memory.verbose", FALSE)
    if (big.memory) 
        for (n in names(ratios)) {
            ratios[[n]] <- matrix.reference(ratios[[n]], backingfile = paste("ratios.", 
                n, sep = ""))
        }
    if (organism == "hsa") 
        rsat.urls[1] <- rsat.urls[2]
    if (!exists("rsat.species") || rsat.species == "?" || is.na(rsat.species)) {
        err <- dlf("data/KEGG/KEGG_taxonomy.txt", "http://baliga.systemsbiology.net/cmonkey/taxonomy.txt")
        if (class(err) != "try-error") {
            tab <- read.delim("data/KEGG/KEGG_taxonomy.txt", 
                sep = "\t", comment = "#", head = F, as.is = T)
            rsat.spec <- as.character(subset(tab, V2 == organism, 
                select = "V4", drop = T))[1]
            rm(tab)
            if (any(strsplit(rsat.spec, "")[[1]] == "(")) 
                rsat.spec <- gsub("\\s\\(.*\\)", "", rsat.spec)
        }
        else {
            rsat.spec <- "?"
        }
        rsat.spec <- gsub(" ", "_", rsat.spec, fixed = T)
        kegg.spec <- rsat.spec
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
            message("Could not find correct organism for RSAT... will try to guess...")
            max.dist <- 0.5
            vals <- rep("", 3)
            while (length(vals) > 1) {
                vals <- agrep(rsat.spec, tmp, ignore = T, max.dist = max.dist, 
                  val = T)
                max.dist <- max.dist - 0.01
                if (length(vals) <= 0) {
                  max.dist <- max.dist + 0.02
                  vals <- agrep(rsat.spec, tmp, ignore = T, max.dist = max.dist, 
                    val = T)
                  break
                }
            }
            if (length(vals) > 1) {
                rsat.spec <- sapply(strsplit(vals, "[<>/]"), 
                  "[", 8)
                message("Found ", length(rsat.spec), " matches...")
                rsat.spec <- rsat.spec[menu(rsat.spec, graphics = F, 
                  title = "Please choose one.")]
            }
            if (length(vals) == 1) {
                rsat.spec <- strsplit(vals, "[<>/]")[[1]][8]
                message("Found one match: ", rsat.spec, " ...")
                message("If this is not correct, you're not quite out of luck -- set the 'rsat.species' parameter manually.")
            }
        }
        set.param("rsat.species", rsat.spec, override = T)
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
        if (cog.o == "") 
            cog.o <- "?"
        set.param("cog.org", cog.o, override = T)
        rm(cog.o, tmp)
    }
    else {
        set.param("cog.org", cog.org)
    }
    message("Organism is ", organism, " ", cog.org, " ", rsat.species, 
        " ", taxon.id)
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
        message("Organism is a eukaryote; presuming there are no operons.")
        set.param("is.eukaryotic", TRUE, override = T)
        set.param("operon.shift", FALSE, override = T)
        set.param("discard.genome", TRUE, override = T)
        set.param("recalc.bg", FALSE, override = T)
        if ("operons" %in% names(net.weights)) {
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
    if (sum(row.weights, na.rm = T) > 0) 
        row.weights <- row.weights/sum(row.weights, na.rm = T)
    if (sum(mot.weights, na.rm = T) > 0) 
        mot.weights <- mot.weights/sum(mot.weights, na.rm = T)
    for (i in c("n.motifs", "meme.cmd", "mast.cmd", "meme.iters", 
        "operon.shift", "bg.order", "motif.upstream.search", 
        "motif.upstream.scan")) {
        v <- get(i)
        if (all(names(mot.weights) %in% names(v))) 
            next
        if (is.vector(v) && length(v) > 1) 
            v <- list(`1` = v)
        names(v) <- names(mot.weights)[1]
        for (n in names(mot.weights)[!names(mot.weights) %in% 
            names(v)]) {
            if (is.list(v)) 
                v[[n]] <- v[[1]]
            else if (is.vector(v)) 
                v[n] <- v[1]
            names(v)[length(v)] <- names(mot.weights)[length(v)]
        }
        assign(i, v)
    }
    rm(v)
    if (!is.null(env)) 
        for (i in ls()) {
            if (i %in% c("i", "env")) 
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
            message("Organism is ", organism, " ", cog.org, " ", 
                rsat.species, " ", taxon.id)
        }
        genome.info$operons <- NULL
        if ((operon.shift || "operons" %in% names(net.weights)) && 
            !no.genome.info) {
            tmp.operons <- try(get.operon.predictions("microbes.online"))
            if (class(tmp.operons) == "try-error") {
                message("Could not fetch operons file. Assuming it doesn't exist (eukaryote?)")
                set.param("is.eukaryotic", TRUE, override = T)
                set.param("operon.shift", FALSE, override = T)
                operon.shift[1:length(operon.shift)] <- FALSE
                if ("operons" %in% names(net.weights)) {
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
        if (exists("ratios") && !is.null(ratios)) 
            tmp <- toupper(attr(ratios, "rnames"))
        else if (exists("genome.info") && !is.null(genome.info$feature.names)) {
            tmp <- toupper(subset(genome.info$feature.names, 
                type == "primary", select = "names", drop = T))
            if (exists("ratios") && !is.null(ratios)) 
                tmp <- tmp[toupper(tmp) %in% toupper(attr(ratios, 
                  "rnames"))]
        }
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
        prefix <- NA
        if (nch > 0) {
            prefix <- names(which.max(table(substr(tmp, 1, nch))))
            message("Assuming gene/probe names have common prefix '", 
                prefix, "'.")
            genome.info$gene.prefix <- prefix
        }
        else {
            message("Could not find a common gene/probe identifier prefix. This only matters if there's no expression matrix.")
            prefix <- genome.info$gene.prefix <- NA
        }
        if (!is.null(env)) 
            assign("genome.info", genome.info, envir = env)
        if (!exists("ratios") || is.null(ratios)) {
            message("WARNING: No ratios matrix -- will generate an 'empty' one with all annotated ORFs for 'probes'.")
            if (!is.na(prefix)) 
                rows <- unique(as.character(subset(genome.info$feature.names, 
                  grepl(paste("^", prefix, sep = ""), names, 
                    ignore = T, perl = T), select = "names", 
                  drop = T)))
            else rows <- unique(as.character(subset(genome.info$feature.names, 
                type == "primary", select = "names", drop = T)))
            ratios <- list(ratios = t(t(rep(NA, length(rows)))))
            rownames(ratios$ratios) <- rows
            attr(ratios, "rnames") <- sort(unique(rows))
            rm(rows)
            attr(ratios, "nrow") <- length(attr(ratios, "rnames"))
            attr(ratios, "ncol") <- 1
            cat("Ratios: ", attr(ratios, "nrow"), "x", 1, "\n")
        }
        rm(nch, prefix, tmp, qqq)
        if (!no.genome.info && length(mot.weights) > 0) {
            genome.info$all.upstream.seqs <- genome.info$bg.list <- list()
            genome.info$bg.fname <- character()
            for (i in names(mot.weights)) {
                cat("Pre-computing all '", i, "' seqs (", paste(motif.upstream.scan[[i]], 
                  collapse = ", "), ")...\n", sep = "")
                genome.info$all.upstream.seqs[[i]] <- get.sequences(attr(ratios, 
                  "rnames"), seq.type = i, distance = motif.upstream.scan[[i]], 
                  filter = F)
                if (!is.null(env)) 
                  assign("genome.info", genome.info, envir = env)
                message(sum(!attr(ratios, "rnames") %in% names(genome.info$all.upstream.seqs[[i]])), 
                  " probes have no '", i, "' sequence.")
                if (!is.na(bg.order[i])) {
                  cat("Pre-computing '", i, "' residue bg distrib (order=", 
                    bg.order[i], ")...\n", sep = "")
                  tmp.seqs <- if (!is.null(genome.info$all.upstream.seqs[[i]])) 
                    genome.info$all.upstream.seqs[[i]]
                  else get.sequences(attr(ratios, "rnames"), 
                    distance = motif.upstream.search[[i]], seq.type = i, 
                    filter = F)
                  genome.info$bg.fname[i] <- my.tempfile("meme.tmp", 
                    suf = ".bg")
                  capture.output(genome.info$bg.list[[i]] <- mkBgFile(tmp.seqs, 
                    order = bg.order[i], bgfname = genome.info$bg.fname[i], 
                    use.rev.comp = grepl("-revcomp", meme.cmd[i])))
                  rm(tmp.seqs)
                }
                else {
                  message("NOT USING a global sequence background distribution!")
                }
                if (!is.null(env)) 
                  assign("genome.info", genome.info, envir = env)
            }
            if (discard.genome) {
                cat("Clearing genome sequences from memory.\n")
                genome.info$genome.seqs <- NULL
                gc()
            }
        }
        networks <- list()
        if (!is.na(net.iters) && any(net.iters %in% 1:n.iter)) {
            if (length(grep("string", names(net.weights))) > 
                0) {
                if ("string" %in% names(net.weights)) {
                  if (exists("string.links")) {
                    string <- string.links
                  }
                  else {
                    cat("Loading STRING network.\n")
                    string <- get.STRING.links(genome.info$org.id$V1[1])
                  }
                  if (!is.null(string) && nrow(string) > 0) {
                    cat("Read in", nrow(string), "STRING edges; weight =", 
                      net.weights["string"], "\n")
                    string$combined_score <- string$combined_score/max(string$combined_score, 
                      na.rm = T) * 1000
                    string$combined_score <- 1000 * exp(string$combined_score/1000)/exp(1)
                    networks[["string"]] <- string
                  }
                  else {
                    warning("Could not load STRING network. Either", 
                      organism, "is not there or your gene names are not standard.")
                  }
                  rm(string)
                }
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if ("operons" %in% names(net.weights) && !is.null(genome.info$operons)) {
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
                if (!is.null(out.sif) && nrow(out.sif) > 0) {
                  out.sif$combined_score <- rep(1000, nrow(out.sif))
                  colnames(out.sif) <- c("protein1", "protein2", 
                    "combined_score")
                  networks[["operons"]] <- out.sif
                }
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
                0 && !is.null(names(net.weights))) {
                for (i in names(net.weights)) {
                  if (i %in% names(networks)) 
                    next
                  if (file.exists(i)) {
                    cat("Loading sif interactions from file:", 
                      i, "; weight =", net.weights[i], "\n")
                    sif <- load.sif.interactions(i)
                  }
                  else if (exists(i) && !is.null(ncol(get(i))) && 
                    ncol(get(i)) >= 2) {
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
                  if (i %in% names(networks)) 
                    next
                  if (file.exists(i)) {
                    cat("Loading groupings from file:", i, "; weight =", 
                      grouping.weights[i], "\n")
                    sif <- load.sif.interactions(i)
                  }
                  else if (exists(i) && !is.null(ncol(get(i))) && 
                    ncol(get(i)) >= 2) {
                    cat("Using groupings from '", i, "' that exists in memory already; weight = ", 
                      grouping.weights[i], "\n", sep = "")
                    sif <- get(i)
                    if (ncol(sif) == 2) 
                      sif <- cbind(sif, combined_score = rep(1, 
                        nrow(sif)))
                  }
                  colnames(sif) <- c("group", "protein", "combined_score")
                  if (sum(unique(as.character(sif$protein)) %in% 
                    attr(ratios, "rnames")) < sum(unique(as.character(sif$group)) %in% 
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
                    if (j%%100 == 0) 
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
            for (n in names(networks)) {
                nn <- networks[[n]]
                if (nrow(nn) <= 0) {
                  message("WARNING: no edges in network", n, 
                    "... skipping.")
                  if (length(grep(n, seed.method[1])) > 0) {
                    message("ALSO, we have to change the row seeding method from", 
                      seed.method, "to 'kmeans'.")
                    seed.method["rows"] <- "kmeans"
                    set.param("seed.method", seed.method, override = T)
                  }
                  next
                }
                nodes <- unique(c(as.character(nn$protein1), 
                  as.character(nn$protein2)))
                cat(nrow(nn), "edges,", length(nodes), "nodes in network", 
                  n, "\n")
                nn <- subset(nn, as.character(protein1) != as.character(protein2))
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
                      1]], protein2 = colnames(tmp)[tmp2[, 2]], 
                      combined_score = tmp[tmp2])
                    rm(tmp, tmp2)
                    nn <- rbind(nn.new, nn)
                    rm(nn.new)
                    nn <- nn[!duplicated(nn[, c("protein1", "protein2")]), 
                      ]
                  }
                  rm(tmp.nn, dupe.nodes)
                }
                if (exists("ratios") && !is.null(ratios) && !any(nodes %in% 
                  attr(ratios, "rnames"))) {
                  if (median(nchar(nodes)) > median(nchar(attr(ratios, 
                    "rnames"))) && any(substr(nodes, 1, median(nchar(attr(ratios, 
                    "rnames")))) %in% attr(ratios, "rnames"))) {
                    nn$protein1 <- substr(as.character(nn$protein1), 
                      1, median(nchar(attr(ratios, "rnames"))))
                    nn$protein2 <- substr(as.character(nn$protein2), 
                      1, median(nchar(attr(ratios, "rnames"))))
                    nodes <- unique(c(as.character(nn$protein1), 
                      as.character(nn$protein2)))
                  }
                  if (!is.null(genome.info$synonyms)) {
                    rr <- attr(ratios, "rnames")[!attr(ratios, 
                      "rnames") %in% nodes]
                    if (length(rr) > 0) {
                      cat("Reconciling network", n, length(rr), 
                        "node names with probe names...\n")
                      syns <- get.synonyms(rr)
                      mc <- get.parallel(length(syns))
                      is.there <- unlist(mc$apply(syns, function(i) any(i %in% 
                        nodes)))
                      syns <- syns[is.there]
                      nnc1 <- as.character(nn$protein1)
                      nnc2 <- as.character(nn$protein2)
                      nnc1.t <- !nnc1 %in% attr(ratios, "rnames")
                      nnc2.t <- !nnc2 %in% attr(ratios, "rnames")
                      mc <- get.parallel(2)
                      tmp <- mc$apply(1:2, function(ii) {
                        for (i in names(syns)) {
                          if (ii == 1) 
                            nnc1[nnc1.t & nnc1 %in% syns[[i]]] <- i
                          else nnc2[nnc2.t & nnc2 %in% syns[[i]]] <- i
                        }
                        if (ii == 1) 
                          return(nnc1)
                        else return(nnc2)
                      })
                      nnc1 <- tmp[[1]]
                      nnc2 <- tmp[[2]]
                      rm(tmp, nnc1.t, nnc2.t)
                      cat(sum(!is.there), "probes have no nodes in", 
                        n, "network (but", sum(attr(ratios, "rnames") %in% 
                          nodes, na.rm = T) + sum(is.there), 
                        "do)\n")
                      nn$protein1 <- nnc1
                      nn$protein2 <- nnc2
                      tmp <- nnc1 %in% attr(ratios, "rnames") & 
                        nnc2 %in% attr(ratios, "rnames")
                      nn <- subset(nn, tmp == TRUE)
                      rm(tmp, syns, is.there, nnc1, nnc2, nnc1.t, 
                        nnc2.t, tmp, rr, i)
                    }
                  }
                }
                else {
                  cat(sum(!attr(ratios, "rnames") %in% nodes), 
                    "probes have no nodes in", n, "network (but", 
                    sum(attr(ratios, "rnames") %in% nodes, na.rm = T), 
                    "do)\n")
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
                ms <- max(sums[sums > 0], na.rm = T)
                if (length(sums) > 0 && !is.na(ms)) 
                  for (n in names(networks)) if (sums[n] > 0) 
                    networks[[n]]$combined_score <- networks[[n]]$combined_score/sums[n] * 
                      ms
                rm(n, sums, ms)
            }
            if (!is.null(env)) 
                assign("networks", networks, envir = env)
            if (!is.null(names(net.weights))) 
                names(net.weights) <- basename(names(net.weights))
        }
        if (big.memory) 
            networks <- list.reference(networks, file = sprintf("%s/networks", 
                cmonkey.filename), type = "RDS")
        if (!no.genome.info && cog.org != "" && cog.org != "?" && 
            !is.null(cog.org)) {
            cat("Loading COG functional codes (for plotting), org. code", 
                cog.org, ": trying NCBI whog file...\n")
            genome.info$cog.code <- get.COG.code(cog.org)
        }
        cat(sum(!is.na(genome.info$cog.code)), "genes have a COG code (", 
            if (is.null(genome.info$cog.code)) 
                attr(ratios, "nrow")
            else sum(is.na(genome.info$cog.code)), "do not)\n")
        if (big.memory) 
            genome.info <- list.reference(genome.info, file = sprintf("%s/genome.info", 
                cmonkey.filename), type = "RDS")
    }
    iter <- 0
    meme.scores <- clusterStack <- list()
    for (i in names(mot.weights)) {
        meme.scores[[i]] <- list()
        meme.scores[[i]][[k.clust + 1]] <- ""
    }
    stats <- row.scores <- col.scores <- mot.scores <- net.scores <- NULL
    if (!exists("favorite.cluster")) 
        favorite.cluster <- function() min(which(sapply(1:k.clust, 
            function(k) length(get.rows(k))) > cluster.rows.allowed[1] * 
            2))
    row.scaling <- extend.vec(row.scaling)
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
        if (i %in% c("i", "env")) 
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
        try(env$cm.func.each.iter(), silent = T)
    }
    cat("INITIALIZATION IS COMPLETE.\n")
    env$iter <- env$iter + 1
    invisible(env)
}
cmonkey.one.iter <-
function (env, dont.update = F, ...) 
{
    env <- env$update.all.clusters(env, dont.update = F, ...)
    row.memb <- sapply(1:k.clust, function(k) attr(ratios, "rnames") %in% 
        get.rows(k))
    if (is.vector(row.memb)) 
        row.memb <- t(row.memb)
    rownames(row.memb) <- attr(ratios, "rnames")
    col.memb <- sapply(1:k.clust, function(k) attr(ratios, "cnames") %in% 
        get.cols(k))
    if (is.vector(col.memb)) 
        col.memb <- t(col.memb)
    rownames(col.memb) <- attr(ratios, "cnames")
    if (iter %in% stats.iters) {
        env$stats <- rbind(env$stats, env$get.stats())
        cat(organism, as.matrix(env$stats[nrow(env$stats), ]), 
            "\n")
    }
    else {
        cat(sprintf("==> %04d %.3f %.3f %.3f\n", iter, mean(env$row.scores[, 
            ][row.memb], na.rm = T), if (!is.null(env$mot.scores)) 
            mean(env$mot.scores[, ][row.memb & env$mot.scores[, 
                ] < 0], na.rm = T, trim = 0.05)
        else NA, if (!is.null(env$net.scores)) 
            mean(env$net.scores[, ][row.memb], na.rm = T, trim = 0.05)
        else NA))
    }
    if (!is.na(plot.iters) && iter %in% plot.iters) {
        env$plotStats(iter, plot.clust = env$favorite.cluster(), 
            new.dev = T)
    }
    if (exists("cm.func.each.iter")) 
        try(cm.func.each.iter(), silent = T)
    if (any(cm.script.each.iter != "")) {
        for (f in cm.script.each.iter) {
            if (file.exists(f) && file.info(f)$size > 1) {
                tmp <- readLines(f)
                if (all(substr(tmp, 1, 1) == "#")) 
                  next
                if (tmp[1] != "## QUIET") 
                  cat("Sourcing the script '", f, "' ...\n", 
                    sep = "")
                try(source(f, echo = tmp[1] != "## QUIET", local = T), 
                  silent = T)
                rm(tmp)
            }
        }
    }
    if (get.parallel()$mc) {
        if (getDoParName() == "doMC") {
            chld <- multicore::children()
            if (length(chld) > 0) {
                try({
                  multicore::kill(chld)
                  tmp <- multicore::collect(chld)
                }, silent = T)
            }
        }
        else if (getDoParName() == "doSNOW" && "data" %in% ls(pos = foreach:::.foreachGlobals)) {
            cl <- get("data", pos = foreach:::.foreachGlobals)
            if (!is.null(data)) 
                stopCluster(cl)
        }
    }
    if (!dont.update) 
        env$iter <- env$iter + 1
    invisible(env)
}
cmonkey.re.seed <-
function (env) 
{
    if (!exists("rnd.seed", envir = env$cmonkey.params)) {
        op <- options(digits.secs = 10)
        tmp.time <- as.character(Sys.time())
        options(op)
        rm(op)
        tmp.rnd.seed <- as.integer(substr(gsub("[-:. ]", "", 
            tmp.time), 12, 20))
        cat("RESETTING RANDOM SEED: ")
        env$set.param("date.run", tmp.time, env$cmonkey.params)
        env$date.run <- env$cmonkey.params$date.run
        env$set.param("rnd.seed", tmp.rnd.seed, env$cmonkey.params)
        env$rnd.seed <- env$cmonkey.params$rnd.seed
        set.seed(env$rnd.seed)
        rm(tmp.rnd.seed)
    }
    if (!is.null(env$ratios) && attr(env$ratios, "ncol") > 1) {
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
    env$clusterStack <- lapply(1:env$k.clust, function(k) list(rows = rownames(which(tmp$rm == 
        k, arr = T)), cols = rownames(which(tmp$cm == k, arr = T))))
    attr(env$clusterStack, "iter") <- env$iter - 1
    invisible(env)
}
col.let <-
c("A", "C", "G", "T")
compare.pssms <-
function (pssm1, pssm2, rev.comp = F, weight = F, min.ov = 6, 
    score = "cor") 
{
    getEntropy = function(pssm) {
        pssm[pssm == 0] <- 1e-05
        entropy <- apply(pssm, 1, function(i) -sum(i * log2(i)))
        return(entropy)
    }
    if (rev.comp) 
        pssm2 <- pssm2[nrow(pssm2):1, 4:1]
    if (score == "cor" || score == "pearson") 
        cors <- cor(t(pssm1), t(pssm2), method = "pearson")
    else if (score == "spearman") 
        cors <- cor(t(pssm1), t(pssm2), method = "spearman")
    else if (score == "ed") 
        cors <- -as.matrix(dist(rbind(pssm1, pssm2), "euclidean"))[1:nrow(pssm1), 
            nrow(pssm1) + 1:nrow(pssm2)]
    if (weight) {
        score1 <- (2 - getEntropy(pssm1))
        score2 <- (2 - getEntropy(pssm2))
        scores <- outer(score1, score2, fun = "min")/2
        cors <- cors * scores
    }
    rev <- FALSE
    if (ncol(cors) > nrow(cors)) {
        rev <- TRUE
        cors <- t(cors)
    }
    min.ind <- ncol(cors)
    max.ind <- nrow(cors)
    vec <- 1:min.ind
    covs <- do.call(rbind, lapply((-max.ind + min.ov):(max.ind - 
        min.ov), function(i) {
        vec2 <- vec + i
        vec2a <- vec2 > 0 & vec2 <= max.ind
        vec2 <- vec2[vec2a]
        if (length(vec2) < min.ov) 
            return(NULL)
        c(i, length(vec2), mean(cors[cbind(vec2, vec[vec2a])]))
    }))
    if (!rev) 
        covs[, 1] <- -covs[, 1]
    return(covs)
}
consolidate.duplicate.clusters <-
function (row.membership, col.membership, scores = r.scores, 
    cor.cutoff = 0.9, n.cutoff = 5, motif = F, seq.type = "upstream meme") 
{
    row.m <- row.membership
    ms <- meme.scores
    cors <- id.duplicate.clusters(scores, cor.cutoff)
    if (nrow(cors) <= 0) 
        return(invisible(list(r = row.m, ms = meme.scores, scores = scores)))
    cr <- max(cors[, 3], na.rm = T)
    n.cut <- 1
    while (cr > cor.cutoff && !is.infinite(cr) && n.cut <= n.cutoff) {
        tmp <- cors[which(cors[, 3] == cr), 1:2]
        if (any(get.rows(tmp[1]) %in% get.rows(tmp[2]))) {
            ev1 <- if (is.null(meme.scores[[seq.type]][[tmp[1]]]$meme.out)) 
                Inf
            else min(sapply(meme.scores[[seq.type]][[tmp[1]]]$meme.out, 
                "[[", "e.value"), na.rm = T)
            ev2 <- if (is.null(meme.scores[[seq.type]][[tmp[2]]]$meme.out)) 
                Inf
            else min(sapply(meme.scores[[seq.type]][[tmp[2]]]$meme.out, 
                "[[", "e.value"), na.rm = T)
            if (!(is.infinite(ev1) && is.infinite(ev2)) && ev2 < 
                ev1) 
                tmp <- tmp[c(2, 1)]
            else if (length(get.rows(tmp[1])) < length(get.rows(tmp[2]))) 
                tmp <- tmp[c(2, 1)]
            row.m[row.m == tmp[2]] <- tmp[1]
            cat("MERGING:", tmp, "\t", length(get.rows(tmp[1])), 
                length(get.rows(tmp[2])), "\t", length(unique(c(get.rows(tmp[1]), 
                  get.rows(tmp[2])))), "\t", cr, "\n")
            scores[, tmp[2]] <- NA
            for (tt in names(mot.weights)) {
                ms[[tt]][[tmp[2]]] <- list(iter = iter)
                if (motif && sum(!get.rows(tmp[1]) %in% get.rows(tmp[2])) > 
                  0) 
                  ms[[tt]][[tmp[1]]] <- try(meme.one.cluster(tmp[1], 
                    verbose = T, consens = meme.consensus, seq.type = tt))
            }
            n.cut <- n.cut + 1
        }
        cors[which(cors[, 3] == cr), ] <- NA
        cr <- max(cors[, 3], na.rm = T)
    }
    invisible(list(r = row.m, ms = ms, scores = scores))
}
cosmo.one.cluster <-
function (k, seq.type = "upstream cosmo", n.motifs = 2, minW = 6, 
    maxW = 24, model = "ZOOPS", verbose = F, unlink = T, ...) 
{
    require(cosmo)
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    seqs <- get.sequences(rows, seq.type = seq.type, ...)
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    bgs <- genome.info$all.upstream.seqs[[seq.type]]
    ss <- list()
    for (s in names(seqs)) ss[[s]] <- list(desc = s, seq = seqs[s])
    bgsl <- list()
    for (s in names(bgs)) bgsl[[s]] <- list(desc = s, seq = bgs[s])
    bgm <- bgModel(bgsl)
    res <- list()
    for (i in 1:n.motifs) {
        res[[i]] <- cosmo(seqs = ss, minW = minW, maxW = maxW, 
            models = model, transMat = bgm$transMat, minSites = length(seqs)/2, 
            maxSites = 2 * length(seqs), wCrit = "eval", intCrit = "eval", 
            silent = !verbose, ...)
        pssm <- t(attr(attr(res[[i]], "pwm"), "pwm")[col.let, 
            ])
        if (i < n.motifs) {
            sites <- attr(res[[i]], "motifs")
            gene <- attr(sites, "seq")
            start <- attr(sites, "pos")
            for (i in 1:length(start)) ss[[gene[i]]]$seq <- substr(ss[[gene[i]]]$seq, 
                start[i], start[i] + nrow(pssm) - 1) <- paste(rep("N", 
                nrow(pssm)), collapse = "")
        }
    }
    m.in <- character()
    for (i in 1:length(res)) {
        pssm <- t(attr(attr(res[[i]], "pwm"), "pwm")[col.let, 
            ])
        m.in <- c(m.in, pssm.motif.lines(pssm, id = sprintf("cosmo_%d", 
            i), header = (i == 1)))
    }
    all.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    mast.out <- runMast(m.in, names(all.seqs), all.seqs, bg.list = genome.info$bg.list[[seq.type]], 
        unlink = T, verbose = verbose)
    pv.ev <- get.pv.ev.single(mast.out, rows)
    meme.out <- list()
    for (ii in 1:length(res)) {
        pssm <- t(attr(attr(res[[ii]], "pwm"), "pwm")[col.let, 
            ])
        pssm <- pssm + max(pssm, na.rm = T)/100
        for (i in 1:nrow(pssm)) pssm[i, ] <- pssm[i, ]/sum(pssm[i, 
            ], na.rm = T)
        sites <- attr(res[[ii]], "motifs")
        posns <- data.frame(gene = attr(sites, "seq"), strand = ifelse(attr(sites, 
            "orient") < 0, "-", "+"), start = attr(sites, "pos"), 
            p.value = 1 - attr(sites, "prob") + 1e-05, site = attr(sites, 
                "motif"))
        meme.out[[ii]] <- list(width = nrow(pssm), sites = nrow(sites), 
            llr = attr(res, "sel")["Width", "critVal"], e.value = attr(res, 
                "sel")["Width", "critVal"], pssm = pssm, posns = posns)
    }
    attr(meme.out, "is.pal") <- FALSE
    invisible(list(k = k, cosmo.out = res, meme.out = meme.out, 
        pv.ev = pv.ev))
}
desymmetrize.tomtom.results <-
function (tt.out) 
{
    mot.names <- unique(c(paste(tt.out$biclust1, tt.out$motif1, 
        sep = "_"), paste(tt.out$biclust2, tt.out$motif2, sep = "_")))
    mot.names.2 <- cbind(paste(tt.out$biclust1, tt.out$motif1, 
        sep = "_"), paste(tt.out$biclust2, tt.out$motif2, sep = "_"))
    mot.names.2a <- cbind(mot.names.2[, 2], mot.names.2[, 1])
    tmp <- matrix(NA, nrow = length(mot.names), ncol = length(mot.names))
    rownames(tmp) <- colnames(tmp) <- mot.names
    lt <- lower.tri(tmp)
    tmp[mot.names.2] <- tt.out$p.value
    tmp2 <- cbind(tmp[lt], t(tmp)[lt])
    tmp2[is.na(tmp2)] <- Inf
    tmp[lt] <- apply(tmp2, 1, min, na.rm = T)
    tmp[is.infinite(tmp)] <- NA
    rm(tmp2)
    tmp.good <- !is.na(tmp[lt])
    tmp2 <- tmp * NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$biclust1
    out <- data.frame(biclust1 = tmp2[lt][tmp.good])
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$motif1
    out <- cbind(out, data.frame(motif1 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$resid1
    out <- cbind(out, data.frame(resid1 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$e.value1
    out <- cbind(out, data.frame(e.value1 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$biclust2
    out <- cbind(out, data.frame(biclust2 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$motif2
    out <- cbind(out, data.frame(motif2 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$resid2
    out <- cbind(out, data.frame(resid2 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$e.value2
    out <- cbind(out, data.frame(e.value2 = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$offset
    out <- cbind(out, data.frame(offset = tmp2[lt][tmp.good]))
    out <- cbind(out, data.frame(p.value = tmp[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$q.value
    out <- cbind(out, data.frame(q.value = tmp2[lt][tmp.good]))
    tmp2[, ] <- NA
    tmp2[mot.names.2] <- tmp2[mot.names.2a] <- tt.out$overlap
    out <- cbind(out, data.frame(overlap = tmp2[lt][tmp.good]))
    rm(tmp2)
    if ("consensus1" %in% colnames(tt.out) || "consensus2" %in% 
        colnames(tt.out) || is.character(tt.out$orientation)) {
        tmp.c <- matrix("", nrow = length(mot.names), ncol = length(mot.names))
        rownames(tmp.c) <- colnames(tmp.c) <- mot.names
        if ("consensus1" %in% colnames(tt.out)) {
            tmp.c[mot.names.2] <- tmp.c[mot.names.2a] <- as.character(tt.out$consensus1)
            out <- cbind(out, data.frame(consensus1 = tmp.c[lt][tmp.good]))
        }
        if ("consensus2" %in% colnames(tt.out)) {
            tmp.c[, ] <- ""
            tmp.c[mot.names.2] <- tmp.c[mot.names.2a] <- as.character(tt.out$consensus2)
            out <- cbind(out, data.frame(consensus2 = tmp.c[lt][tmp.good]))
        }
        if (is.factor(tt.out$orientation)) 
            tt.out$orientation <- as.character(tt.out$orientation)
        if (is.character(tt.out$orientation)) {
            tmp.c[, ] <- ""
            tmp.c[mot.names.2] <- tmp.c[mot.names.2a] <- as.character(tt.out$orientation)
            out <- cbind(out, data.frame(orientation = tmp.c[lt][tmp.good]))
        }
        rm(tmp.c)
    }
    out <- subset(out, !is.na(p.value))
    out <- out[order(out$p.value, out$q.value), ]
    out
}
dlf <-
function (f, url, msg = NULL, mode = "wb", quiet = F, ...) 
{
    err <- 0
    if (mode == "ab" || !file.exists(f) || file.info(f)$size == 
        0) {
        if (!file.exists(dirname(f))) 
            try(dir.create(dirname(f), recursive = T))
        if (!is.null(msg)) 
            cat(msg, "\n")
        err <- try(download.file(url, destfile = f, mode = mode, 
            quiet = quiet, ...))
    }
    closeAllConnections()
    err
}
extend.vec <-
function (v, n = n.iter) 
{
    if (length(v) < n) 
        v <- c(v, rep(v[length(v)], n.iter - length(v)))
    v
}
ffify.env <-
function (env) 
{
    for (i in c("row.scores", "mot.scores", "net.scores")) {
        if (exists(i, envir = env)) {
            tmp <- matrix.reference(env[[i]], backingfile = i, 
                backingpath = env$cmonkey.filename)
            assign(i, tmp, envir = env)
        }
    }
    for (i in names(env$ratios)) {
        env$ratios[[i]] <- matrix.reference(env$ratios[[i]], 
            backingfile = paste("ratios.", i, sep = ""), backingpath = env$cmonkey.filename)
    }
    for (i in names(env$meme.scores)) {
        file <- paste(env$cmonkey.filename, "/meme.scores.", 
            i, sep = "")
        if (exists("meme.scores", envir = env)) {
            env$meme.scores[[i]] <- list.reference(env$meme.scores[[i]], 
                file)
        }
    }
    for (i in c("clusterStack", "genome.info", "networks")) {
        file <- paste(env$cmonkey.filename, "/", i, sep = "")
        if (exists(i, envir = env)) 
            env[[i]] <- list.reference(env[[i]], file)
    }
    invisible(env)
}
filter.sequences <-
function (seqs, start.stops = NULL, seq.type = paste(c("upstream", 
    "upstream.noncod", "upstream.noncod.same.strand", "downstream", 
    "gene")[1], "meme"), distance = motif.upstream.search[[seq.type]], 
    uniquify = T, remove.repeats = T, remove.atgs = T, mask.overlapping.rgns = F, 
    blast.overlapping.rgns = F, verbose = F, ...) 
{
    if (uniquify) 
        seqs <- seqs[!get.dup.seqs(seqs)]
    if (remove.repeats && length(grep("NNNNNN", seqs)) <= 1) {
        if (verbose) 
            cat("Removing low-complexity regions from sequences.\n")
        seqs.new <- remove.low.complexity(seqs, seq.type = seq.type)
        if (length(seqs.new) == length(seqs)) 
            seqs <- seqs.new
        else warning("Remove low complexity failed - skipping!")
        rm(seqs.new)
    }
    if (remove.atgs && any(distance < 0)) {
        tmp <- names(seqs)
        substr(seqs, distance[2] + 1, distance[2] + 4) <- "NNNN"
        names(seqs) <- tmp
    }
    if (mask.overlapping.rgns) {
        if (is.null(start.stops)) 
            start.stops <- attr(seqs, "start.stops")
        if (!is.null(start.stops)) {
            overlaps <- apply(start.stops, 1, function(i) subset(start.stops, 
                i[4] == contig & (i[1] >= start & i[1] <= end) | 
                  (i[2] >= start & i[2] <= end)))
            overlaps <- lapply(names(overlaps), function(g) overlaps[[g]][rownames(overlaps[[g]]) != 
                g, ])
            names(overlaps) <- rownames(start.stops)
            is.overlapping <- sapply(overlaps, nrow)
            overlaps <- overlaps[is.overlapping > 0]
            for (i in names(overlaps)) {
                if (nrow(overlaps[[i]]) <= 0) 
                  next
                seq1 <- seqs[i]
                if (start.stops[i, 3] == "R") 
                  seq1 <- rev.comp(seq1)
                ss1 <- sapply(20:nchar(seq1), function(i) substr(seq1, 
                  1, i))
                ss2 <- sapply(1:(nchar(seq1) - 20), function(i) substr(seq1, 
                  i, nchar(seq1)))
                for (j in 1:nrow(overlaps[[i]])) {
                  seq2 <- seqs[rownames(overlaps[[i]])[j]]
                  if (overlaps[[i]][j, 3] == "R") 
                    seq2 <- rev.comp(seq2)
                  g1 <- sapply(sapply(ss1, grep, seq2), length)
                  rgn <- c(1, nchar(seq1))
                  if (all(g1 > 0)) {
                  }
                  else if (any(g1 > 0)) {
                    ind <- which(diff(g1) != 0)
                    rgn <- c(1, ind - 1)
                  }
                  else {
                    g2 <- sapply(sapply(ss2, grep, seq2), length)
                    if (any(g2 > 0)) {
                      ind <- which(diff(g2) != 0)
                      rgn <- c(ind + 1, nchar(seq1))
                    }
                  }
                  if (verbose) 
                    cat(sprintf("Masking region %d-%d of sequence %s (%s)\n", 
                      rgn[1], rgn[2], i, rownames(overlaps[[i]])[j]))
                  substr(seq1, rgn[1], rgn[2]) <- paste(rep("N", 
                    rgn[2] - rgn[1] + 1), collapse = "")
                  seq <- seq1
                  if (start.stops[i, 3] == "R") 
                    seq <- rev.comp(seq)
                  seqs[i] <- seq
                  other.ov <- rownames(overlaps[[i]])[j]
                  overlaps[[other.ov]] <- overlaps[[other.ov]][rownames(overlaps[[other.ov]]) != 
                    i, , drop = F]
                }
            }
        }
    }
    if (blast.overlapping.rgns && exists("blast.match.seqs") && 
        file.exists(sprintf("%s/blastall", progs.dir))) {
        out <- blast.match.seqs(seqs)
        while (nrow(out) > 0) {
            out <- subset(out, alignment.length > 30 & X..identity > 
                80)
            if (nrow(out) <= 0) 
                break
            for (i in 1:nrow(out)) {
                seq1 <- strsplit(seqs[as.character(out$Query.id[i])], 
                  "")[[1]]
                seq2 <- strsplit(seqs[as.character(out$Subject.id[i])], 
                  "")[[1]]
                st.st.1 <- c(out$q..start[i], out$q..end[i])
                st.st.2 <- c(out$s..start[i], out$s..end[i])
                n.n.1 <- sum(seq1[st.st.1[1]:st.st.1[2]] == "N")
                if (n.n.1 >= max(st.st.1) - min(st.st.1) + 1) 
                  next
                n.n.2 <- sum(seq2[st.st.2[1]:st.st.2[2]] == "N")
                if (n.n.2 >= max(st.st.2) - min(st.st.2) + 1) 
                  next
                if (n.n.1 == 0 && n.n.2 == 0) {
                  if (verbose) 
                    cat(sprintf("Masking (BLAST) region %d-%d of sequence %s (%s)\n", 
                      st.st.2[1], st.st.2[2], out$Subject.id[i], 
                      out$Query.id[i]))
                  seq2[st.st.2[1]:st.st.2[2]] <- "N"
                  seqs[as.character(out$Subject.id[i])] <- paste(seq2, 
                    collapse = "")
                }
                else if (n.n.1 > n.n.2) {
                  if (verbose) 
                    cat(sprintf("Masking region %d-%d of sequence %s (%s)\n", 
                      st.st.1[1], st.st.1[2], out$Query.id[i], 
                      out$Target.id[i]))
                  seq1[st.st.1[1]:st.st.1[2]] <- "N"
                  seqs[as.character(out$Query.id[i])] <- paste(seq1, 
                    collapse = "")
                }
                else if (n.n.2 >= n.n.1) {
                  if (verbose) 
                    cat(sprintf("Masking (BLAST) region %d-%d of sequence %s (%s)\n", 
                      st.st.2[1], st.st.2[2], out$Subject.id[i], 
                      out$Query.id[i]))
                  seq2[st.st.2[1]:st.st.2[2]] <- "N"
                  seqs[as.character(out$Subject.id[i])] <- paste(seq2, 
                    collapse = "")
                }
            }
            out <- blast.match.seqs(seqs)
        }
    }
    if (!is.null(start.stops)) 
        attr(seqs, "start.stops") <- start.stops[names(seqs), 
            , drop = F]
    seqs
}
filter.updated.memberships <-
function (row.membership, col.membership, rr.scores, cc.scores, 
    quant.cutoff = c(rows = 0, cols = 0)) 
{
    rm <- row.membership
    if (quant.cutoff["rows"] > 0) {
        qc <- quantile(rr.scores[, ][row.memb[, ] == 1], prob = quant.cutoff["rows"])
        for (i in 1:nrow(rm)) {
            tmp <- which(rm[i, ] != 0)
            rm[i, tmp[rr.scores[i, rm[i, tmp]] < qc]] <- 0
        }
    }
    cm <- col.membership
    if (quant.cutoff["cols"] > 0) {
        qc <- quantile(cc.scores[, ][col.memb[, ] == 1], prob = quant.cutoff["cols"])
        for (i in 1:nrow(cm)) {
            tmp <- which(cm[i, ] != 0)
            cm[i, tmp[cc.scores[i, cm[i, tmp]] < qc]] <- 0
        }
    }
    NULL
}
foreach.register.backend <-
function (par) 
{
    if (!require(foreach)) 
        return(NULL)
    if (par > 1 && require(doMC, quietly = T)) 
        registerDoMC(cores = par)
    else registerDoSEQ()
}
gadem.one.cluster <-
function (k, seq.type = names(mot.weights)[1], n.motifs = 2, 
    verbose = F, unlink = T, ...) 
{
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    seqs <- get.sequences(rows, seq.type = seq.type, ...)
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    bgs <- genome.info$all.upstream.seqs[[seq.type]]
    fname <- my.tempfile("gadem_")
    cat(paste(">", names(seqs), "\n", seqs, sep = ""), file = fname, 
        sep = "\n")
    bgtab <- data.frame(names(genome.info$bg.list[[seq.type]]), 
        sapply(genome.info$bg.list[[seq.type]], "[", 1))
    rownames(bgtab) <- NULL
    colnames(bgtab) <- c("V1", "V2")
    bgtab <- bgtab[-1, ]
    bgtab$V1 <- tolower(bgtab$V1)
    bgfname <- my.tempfile("gadem_bg_")
    write.table(bgtab, bgfname, sep = "\t", quote = F, col.names = F, 
        row.names = F)
    outfname <- my.tempfile("gadem_out_")
    cmd <- sprintf("%s/gadem -fseq %s -minN %d -bOrder %d -fbm %s -fout %s -pgf 0 -verbose %d", 
        progs.dir, fname, floor(sqrt(length(seqs))), bg.order[seq.type], 
        bgfname, outfname, if (verbose) 
            1
        else 0)
    paste(cmd, "-pgf 0 -gen 20 -pop 100 -em 200 -fEM 0.9 -maxgap 18 -useScore 1")
    if (verbose) 
        print(cmd)
    out <- system(cmd, intern = T, ignore.stderr = !verbose)
    out <- readLines(outfname)
    if (unlink) 
        unlink(c(fname, bgfname, outfname))
    out
}
get.COG.code <-
function (org, rows = attr(ratios, "rnames")) 
{
    up.rows <- toupper(rows)
    out <- rep("-", length(rows))
    names(out) <- up.rows
    fname <- "data/COG_whog.txt"
    err <- dlf(fname, "ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog", 
        "Fetching COG codes from NCBI...")
    lines <- readLines(gzfile(fname))
    closeAllConnections()
    hits <- grep(paste(org, "\\:|COG", sep = ""), lines)
    hpy.hits <- grep(paste(org, "\\:", sep = ""), lines[hits])
    if (length(hpy.hits) <= 0) 
        return(NULL)
    genes <- gsub(paste("\\s+", org, "\\:\\s+", sep = ""), "", 
        lines[hits][hpy.hits], perl = T)
    cogs <- lines[hits][hpy.hits - 1]
    cog.codes <- sapply(strsplit(cogs, "[\\s+\\[\\]]", perl = T), 
        "[", 2)
    cog.codes <- substr(cog.codes, 1, 1)
    genes <- toupper(genes)
    mc <- get.parallel(length(genes))
    tmp <- mc$apply(1:length(genes), function(i) {
        gn <- strsplit(genes[i], " ")[[1]]
        if (length(gn) <= 0) 
            next
        gn <- gn[!is.na(gn)]
        if (!all(gn %in% up.rows)) 
            gn <- toupper(unlist(get.synonyms(gn, ignore = T)))
        if (sum(gn %in% up.rows) <= 0) 
            return(character())
        gn <- gn[gn %in% up.rows]
        out[up.rows %in% gn] <- cog.codes[i]
        out
    })
    for (t in tmp) if (length(t) > 0) 
        out[t != "-"] <- t[t != "-"]
    out[out == "-"] <- NA
    names(out) <- rows
    closeAllConnections()
    out
}
get.STRING.links <-
function (org.id = genome.info$org.id$V1[1], all.genes = attr(ratios, 
    "rnames"), score = "score", min.score = 2, string.url = "http://string-db.org/") 
{
    if (file.exists(paste("data/", rsat.species, "/string_links_FALSE_", 
        org.id, ".tab", sep = ""))) {
        string.links <- read.delim(paste("data/", rsat.species, 
            "/string_links_FALSE_", org.id, ".tab", sep = ""), 
            head = T, sep = " ")
        string.links$protein1 <- gsub(paste(org.id, ".", sep = ""), 
            "", string.links$protein1)
        string.links$protein2 <- gsub(paste(org.id, ".", sep = ""), 
            "", string.links$protein2)
        return(string.links)
    }
    file <- sprintf("data/%s/string_links_%s.tab", rsat.species, 
        org.id)
    proc.string.df <- function(file) {
        err <- try(tmp <- unique(read.delim(file, head = F, sep = "")))
        if ("try-catch" %in% class(err)) 
            return(NULL)
        if (!exists("tmp") || nrow(tmp) <= 0) 
            return(NULL)
        tmp2 <- strsplit(as.character(tmp$V15), "[:|]", perl = T)
        tmp2a <- sapply(tmp2, function(i) which(i == score))
        tmp2b <- sapply(1:length(tmp2), function(i) if (length(tmp2a[[i]]) == 
            0) 
            NA
        else as.numeric(tmp2[[i]][tmp2a[[i]] + 1]))
        string.links <- data.frame(protein1 = gsub(paste("string:", 
            org.id, ".", sep = ""), "", tmp$V1), protein2 = gsub(paste("string:", 
            org.id, ".", sep = ""), "", tmp$V2), combined_score = tmp2b)
        string.links <- unique(subset(string.links, !is.na(combined_score)))
        string.links
    }
    string.links <- NULL
    tried <- character()
    if (file.exists(file)) 
        string.links <- proc.string.df(file)
    if (file.exists(sprintf("%s.tried", file))) 
        tried <- readLines(sprintf("%s.tried", file))
    tmp2 <- all.genes %in% tried
    if (!file.exists(file) || (!is.null(string.links) && any(!tmp2))) {
        if (!is.null(string.links)) 
            all.genes <- all.genes[!tmp2]
        id.file <- tempfile()
        options(timeout = 300)
        for (i in seq(1, length(all.genes), by = 100)) {
            ids <- paste(org.id, all.genes[i:min(i + 99, length(all.genes))], 
                sep = ".")
            cat(i, "of", length(all.genes), "\n")
            if (org.id == 3702) {
                ids <- all.genes[i:min(i + 99, length(all.genes))]
                url <- paste(string.url, "api/tsv/resolveList?caller_identity=cMonkey&identifiers=", 
                  URLencode(paste(ids, collapse = "\r"), reserved = T), 
                  sep = "")
                dlf(id.file, url, mode = "wb", quiet = T)
                ids <- unique(as.character(read.delim(id.file)$stringId))
                unlink(id.file)
            }
            url <- paste(string.url, "api/psi-mi-tab/interactionsList?required_score=", 
                min.score, "&caller_identity=cMonkey&network_graph=2&limit=99999&identifiers=", 
                URLencode(paste(ids, collapse = "\r"), reserved = T), 
                sep = "")
            if (!file.exists(file)) 
                dlf(file, url, mode = "wb", msg = "Fetching STRING protein links (piecewise)... this may take a while...", 
                  quiet = T)
            else dlf(file, url, mode = "ab", quiet = T)
        }
        string.links <- proc.string.df(file)
        writeLines(unique(c(all.genes, tried)), sprintf("%s.tried", 
            file))
        options(timeout = 60)
    }
    invisible(string.links)
}
get.STRING.links.OLD <-
function (org.id = genome.info$org.id$V1[1], detailed = T) 
{
    url <- string.links.url
    fname <- strsplit(url, "/")[[1]]
    fname <- sprintf("data/STRING/%s", fname[length(fname)])
    small.fname <- paste("data/", rsat.species, "/string_links_", 
        detailed, "_", org.id, ".tab", sep = "")
    if ((!file.exists(small.fname) || file.info(small.fname)$size <= 
        0)) {
        if (!file.exists(fname)) {
            err <- dlf(fname, url, paste("Fetching STRING protein links file", 
                url, "\nThis will take a while...\n"))
            if (class(err) == "try-error" || !file.exists(fname) || 
                file.info(fname)$size < 1e+09) 
                stop(paste("Whoops, the file was not completely downloaded. Please try to download it yourself from", 
                  string.links.url, "and place it in data/STRING/, then restart cMonkey.\n"))
        }
        cat("Loading organism-specific EMBL STRING interaction links (requires UNIX programs \"gunzip\" and \"grep\")", 
            "...\nUsing local file", fname, "->", small.fname, 
            "\n")
        system.time.limit(paste("gunzip -c ", fname, " | grep -E \"combined_score|^", 
            org.id, ".\" > ", small.fname, sep = ""))
    }
    if (file.exists(small.fname) && file.info(small.fname)$size == 
        0) 
        system.time.limit(paste("gunzip -c ", fname, " | grep -E \"combined_score|^", 
            org.id, ".\" > ", small.fname, sep = ""))
    if (file.exists(small.fname) && file.info(small.fname)$size > 
        0) {
        cat("Loading EMBL STRING interaction links from local file", 
            small.fname, "\n")
        string.links <- read.delim(gzfile(small.fname), sep = " ", 
            head = T)
        string.links$protein1 <- gsub(paste(org.id, ".", sep = ""), 
            "", string.links$protein1)
        string.links$protein2 <- gsub(paste(org.id, ".", sep = ""), 
            "", string.links$protein2)
    }
    url <- string.links.url
    fname <- strsplit(url, "/")[[1]]
    fname <- sprintf("data/STRING/%s", fname[length(fname)])
    dlf(gsub(".gz", "", gsub("protein.links", "species", fname)), 
        gsub(".gz", "", gsub("protein.links", "species", url)))
    closeAllConnections()
    invisible(string.links)
}
get.all.scores <-
structure(function (k, return.scores = F, densities = F, verbose = F, 
    force.motif = F, allow.motif = T, members = F, remove.nas = T, 
    plot = F, ...) 
{
    if (verbose) 
        cat("Computing scores for cluster:", k, "\n")
    row.memb <- attr(ratios, "rnames") %in% get.rows(k)
    names(row.memb) <- attr(ratios, "rnames")
    col.memb <- attr(ratios, "cnames") %in% get.cols(k)
    names(col.memb) <- attr(ratios, "cnames")
    x <- NULL
    for (i in names(ratios)) {
        if (row.weights[i] == 0 || is.na(row.weights[i])) 
            next
        r <- get.row.scores(k, ratios = ratios[[i]], method = row.score.func, 
            ...)
        if (remove.nas) {
            tmp <- is.infinite(r) | is.na(r)
            if (any(tmp)) 
                r[tmp] <- quantile(r[row.memb[names(r)] & !tmp], 
                  0.95)
        }
        x <- cbind(x, r)
        colnames(x)[ncol(x)] <- paste("ratios", i, sep = ".")
    }
    y <- NULL
    for (i in names(ratios)) {
        if (row.weights[i] == 0 || is.na(row.weights[i])) 
            next
        c <- get.col.scores(k, ratios = ratios[[i]], method = col.score.func, 
            ...)
        if (all(is.na(c))) 
            names(c) <- colnames(ratios[[i]])
        if (remove.nas) {
            tmp <- is.infinite(c) | is.na(c)
            if (any(tmp)) 
                c[tmp] <- quantile(c[col.memb[names(c)] & !tmp], 
                  0.95)
        }
        y <- cbind(y, c)
        colnames(y)[ncol(y)] <- paste("ratios", i, sep = ".")
    }
    out.ms <- list()
    for (i in names(mot.weights)) {
        if (mot.weights[i] == 0 || is.na(mot.weights[i])) 
            next
        out <- meme.scores[[i]][[k]]
        if (allow.motif == TRUE && (force.motif == TRUE || force.motif == 
            "run.meme" || (mot.scaling[iter] > 0 && !is.na(meme.iters[[i]][1]) && 
            iter %in% meme.iters[[i]] && exists("genome.info") && 
            !no.genome.info))) {
            out <- motif.one.cluster(k, seq.type = i, verbose = verbose, 
                force = force.motif, ...)
            if (is.null(out) || class(out) == "try-error" || 
                out$k != k || (!is.null(out$iter) && out$iter != 
                iter)) {
                out <- try(motif.one.cluster(k, seq.type = i, 
                  verbose = verbose, force = force.motif, ...))
                if (class(out) == "try-error") 
                  out <- try(motif.one.cluster(k, seq.type = i, 
                    verbose = verbose, force = force.motif, ...))
                if (class(out) == "try-error" || is.null(out) || 
                  out$k != k) {
                  message("ERROR on cluster ", k)
                  out <- list()
                }
                else if (verbose) {
                  cat(iter, k, length(get.rows(k)), seq.type, 
                    "\t")
                }
                if (verbose) {
                  if (is.null(out) || is.null(out$meme.out)) 
                    cat("Inf \n")
                  else {
                    ind <- 1
                    if (!is.null(out$pv.ev)) {
                      if ("p.value" %in% colnames(out$pv.ev[[1]])) 
                        mn <- mean(log10(out$pv.ev[[ind]][rownames(out$pv.ev[[ind]]) %in% 
                          get.rows(k), "p.value"]), na.rm = T)
                      else mn <- mean(log10(out$pv.ev[[ind]]$pvals), 
                        na.rm = T)
                    }
                    else {
                      mn <- "Inf"
                    }
                    cat(k, if (attr(out$meme.out, "is.pal")) 
                      "pal"
                    else "non", sapply(out$meme.out[1:min(3, 
                      length(out$meme.out))], "[[", "e.value"), 
                      mn, "\t", pssm.to.string(out$meme.out[[1]]$pssm), 
                      "\n")
                  }
                }
            }
            out$iter <- iter
            out$k <- k
        }
        m <- get.motif.scores(k, out, ...)
        if (remove.nas) 
            m[is.infinite(m) | is.na(m)] <- 0
        if (!is.null(out) && (is.null(out$iter) || out$iter == 
            iter)) 
            out.ms[[i]] <- out
        x <- cbind(x, m)
        colnames(x)[ncol(x)] <- paste("motif", i, sep = ".")
    }
    for (i in names(networks)) {
        if (net.weights[i] == 0 || is.na(net.weights[i])) 
            next
        if (nrow(subset(networks[[i]], protein1 %in% attr(ratios, 
            "rnames") & protein2 %in% attr(ratios, "rnames"))) <= 
            0) 
            next
        n <- get.network.scores(k, net = networks[[i]])
        if (remove.nas) 
            n[is.infinite(n) | is.na(n)] <- 0
        x <- cbind(x, n)
        colnames(x)[ncol(x)] <- paste("network", i, sep = ".")
    }
    rm(r, c, m, n)
    x.d <- y.d <- NULL
    if (densities || members) {
        p <- rep(1/nrow(x), nrow(x))
        x.d <- apply(x, 2, function(xx) {
            if (!all(is.na(xx)) && !all(xx == xx[1])) {
                fun <- ecdf(c(xx[row.memb], max(xx[row.memb]) + 
                  1e-05))
                p <- 1 - fun(xx)
            }
            p
        })
        rownames(x.d) <- rownames(x)
        p <- rep(1/nrow(y), nrow(y))
        y.d <- apply(y, 2, function(yy) {
            if (!all(is.na(yy)) && !all(yy == yy[1])) {
                fun <- ecdf(c(yy[col.memb], max(yy[col.memb]) + 
                  1e-05))
                p <- 1 - fun(yy)
            }
            p
        })
        rownames(y.d) <- rownames(y)
    }
    scores <- list(r = x, c = y, r.d = x.d, c.d = y.d, ms = out.ms, 
        k = k)
    if (plot) 
        plot.cluster.scores(scores)
    if (members) 
        scores$members <- get.cluster.members(scores, ...)
    if (!return.scores) 
        scores$r <- scores$c <- scores$r.d <- scores$c.d <- NULL
    return(scores)
}, version = 2)
get.biclusters <-
function (genes, conditions, motifs, tfs) 
{
    if (!missing(genes)) {
        out <- lapply(genes, function(g) {
            out <- which(sapply(e$clusterStack, function(i) g %in% 
                i$rows))
            paste("BIC", out, sep = "_")
        })
        names(out) <- genes
    }
    if (!missing(conditions)) {
        out <- lapply(conditions, function(cc) {
            out <- which(sapply(e$clusterStack, function(i) cc %in% 
                i$cols))
            paste("BIC", out, sep = "_")
        })
        names(out) <- conditions
    }
    if (!missing(motifs)) {
        out <- lapply(motifs, function(m) paste("BIC", sapply(strsplit(m, 
            "_"), "[", 2), sep = "_"))
        names(out) <- motifs
    }
    if (!missing(tfs)) {
        out <- lapply(tfs, function(tf) lapply(e.coeffs, function(k) {
            tmp <- k[[1]]$coeffs
            tmp[grep(tf, names(tmp))]
        }))
        for (i in 1:length(out)) {
            out[[i]] <- out[[i]][sapply(out[[i]], length) > 0]
            names(out[[i]]) <- paste("BIC", as.integer(names(out[[i]])), 
                sep = "_")
        }
        names(out) <- tfs
    }
    out
}
get.clust <-
function (k, fill = T, fill.motif = T, seq.type = names(mot.weights), 
    varNorm = F, ...) 
{
    gen.clust <- function(rowNames, colNames = NA, fill = F, 
        motif = F, n.motifs = 3, ...) {
        rowNames <- rowNames[rowNames %in% attr(ratios, "rnames")]
        if (!is.null(colNames) && length(colNames) > 1 && !is.na(colNames)) 
            colNames <- colNames[colNames %in% attr(ratios, "cnames")]
        c.tmp <- list(nrows = length(rowNames), ncols = length(colNames), 
            rows = rowNames, cols = colNames, k = 999, p.clust = 1, 
            e.val = rep(999, n.motifs), resid = {
                out = rep(NA, length(row.weights))
                names(out) <- names(row.weights)
                resid = out
            })
        if (fill && c.tmp$nrows > 0 && c.tmp$ncols > 0 && !all(is.na(colNames))) 
            c.tmp$resid <- cluster.resid(k, names(row.weights), 
                varNorm = varNorm, ...)
        return(c.tmp)
    }
    cols <- get.cols(k)
    rows <- get.rows(k)
    if (length(cols) <= 0) 
        cols <- NA
    clust <- gen.clust(rows, cols, fill = fill, motif = F, n.motifs = max(unlist(n.motifs)))
    clust$k <- k
    if (fill.motif) {
        tmp <- cluster.pclust(k, seq.type)
        clust$e.val <- tmp$e.vals
        clust$p.clust <- tmp$p.clusts
    }
    clust
}
get.cluster.matrix <-
function (rows = NULL, cols = NULL, matrices = names(ratios)) 
{
    if (is.null(rows)) 
        rows <- attr(ratios, "rnames")
    if (is.null(cols)) 
        cols <- attr(ratios, "cnames")
    cols.b <- attr(ratios, "cnames")[attr(ratios, "cnames") %in% 
        cols]
    rats <- matrix(NA, nrow = length(rows), ncol = length(cols.b))
    rownames(rats) <- rows
    colnames(rats) <- cols.b
    cnames <- character()
    for (n in matrices) {
        r.tmp <- ratios[[n]][rows[rows %in% rownames(ratios[[n]])], 
            cols.b[cols.b %in% colnames(ratios[[n]])], drop = F]
        if (is.null(r.tmp) || all(is.na(r.tmp))) 
            next
        if (is.vector(r.tmp)) {
            r.tmp <- t(r.tmp)
            rownames(r.tmp) <- rows
        }
        cnames <- c(cnames, colnames(r.tmp))
        rats[rownames(r.tmp), colnames(r.tmp)] <- r.tmp
        rm(r.tmp)
    }
    rats[, colnames(rats) %in% cnames, drop = F]
}
get.cluster.members <-
function (scores, weights = "calc", pseudo = 0.01, TEMP = "calc", 
    max.change = c(rows = 5, cols = 10), count.power = c(rows = 8, 
        cols = 2)) 
{
    if (TEMP == "calc") 
        TEMP <- seq(0.15, 0.05, length = n.iter)[iter]
    wts <- rep(0, ncol(scores$r.d))
    if (!is.na(weights) && weights == "calc") {
        cn <- sapply(strsplit(colnames(scores$r.d), ".", fixed = T), 
            "[", 1)
        cn2 <- sapply(lapply(strsplit(colnames(scores$r.d), ".", 
            fixed = T), "[", -1), paste, collapse = ".")
        wts[cn == "ratios"] <- row.scaling[iter] * row.weights[cn2[cn == 
            "ratios"]]
        wts[cn == "motif"] <- mot.scaling[iter] * mot.weights[cn2[cn == 
            "motif"]]
        wts[cn == "network"] <- net.scaling[iter] * net.weights[cn2[cn == 
            "network"]]
        wts[is.na(wts)] <- 0
        wts <- wts/sum(wts, na.rm = T)
    }
    if (all(wts == 0 | is.na(wts))) 
        wts[] <- 1
    probs <- apply(scores$r.d, 1, weighted.mean, w = wts, na.rm = T)
    probs <- probs/max(probs, na.rm = T)
    rows <- get.rows(scores$k)
    if (length(rows) > 0 && !all(is.na(rows)) && !is.na(count.power) && 
        !is.na(count.power["rows"])) {
        counts <- table(unlist(lapply(clusterStack, "[[", "rows")))[names(probs)]
        counts[is.na(counts)] <- 0
        names(counts) <- names(probs)
        count.prob <- ppois(counts, n.clust.per.row, lower = F)/ppois(2, 
            n.clust.per.row, lower = F)
        count.prob[counts <= 1] <- 1.1
        count.prob[rows] <- ppois(counts[rows], n.clust.per.row, 
            lower = T)/ppois(n.clust.per.row, n.clust.per.row, 
            lower = T)
        count.prob <- count.prob^count.power["rows"]
    }
    if (FALSE) {
        rows <- unique(c(rows, names(which(probs >= 0.5))))
    }
    else {
        probs.add <- exp(-(1 - probs)/TEMP)
        probs.drop <- exp(-probs/TEMP)
        row.memb <- attr(ratios, "rnames") %in% rows
        names(row.memb) <- attr(ratios, "rnames")
        sample.probs <- ifelse(row.memb, probs.drop, probs.add)
        if (exists("count.prob")) 
            sample.probs <- sample.probs * count.prob
        n.row.per.clust <- n.clust.per.row/k.clust * length(row.memb)
        balance <- sum(sample.probs[row.memb])/sum(sample.probs[!row.memb]) * 
            n.row.per.clust/sum(row.memb)
        if (is.infinite(balance) || is.na(balance)) 
            balance <- 1
        if (balance >= 1) 
            sample.probs[!row.memb] <- sample.probs[!row.memb] * 
                balance
        else sample.probs[row.memb] <- sample.probs[row.memb]/balance
        allowed.moves <- sample.probs[sample.probs > runif(length(sample.probs))]
        if (length(allowed.moves) > max.change["rows"]) 
            allowed.moves <- sample(allowed.moves, max.change["rows"], 
                prob = sample.probs[names(allowed.moves)])
        rows <- c(rows, names(row.memb[names(allowed.moves)][row.memb[names(allowed.moves)] == 
            FALSE]))
        rows <- rows[!rows %in% names(row.memb[names(allowed.moves)][row.memb[names(allowed.moves)] == 
            TRUE])]
    }
    wts <- rep(0, ncol(scores$c.d))
    if (!is.na(weights) && weights == "calc") {
        cn <- sapply(strsplit(colnames(scores$c.d), ".", fixed = T), 
            "[", 1)
        cn2 <- sapply(lapply(strsplit(colnames(scores$c.d), ".", 
            fixed = T), "[", -1), paste, collapse = ".")
        wts[cn == "ratios"] <- row.scaling[iter] * row.weights[cn2[cn == 
            "ratios"]]
        wts[cn == "motif"] <- mot.scaling[iter] * mot.weights[cn2[cn == 
            "motif"]]
        wts[cn == "network"] <- net.scaling[iter] * net.weights[cn2[cn == 
            "network"]]
        wts[is.na(wts)] <- 0
        wts <- wts/sum(wts, na.rm = T)
    }
    if (all(wts == 0 | is.na(wts))) 
        wts[] <- 1
    probs <- apply(scores$c.d, 1, weighted.mean, w = wts, na.rm = T)
    probs <- probs/max(probs, na.rm = T)
    cols <- get.cols(scores$k)
    if (length(cols) > 0 && !all(is.na(cols)) && !is.na(count.power) && 
        !is.na(count.power["cols"])) {
        counts <- table(unlist(lapply(clusterStack, "[[", "cols")))[names(probs)]
        counts[is.na(counts)] <- 0
        names(counts) <- names(probs)
        count.prob <- ppois(counts, n.clust.per.col, lower = F)/ppois(2, 
            n.clust.per.col, lower = F)
        count.prob[counts <= 1] <- 1.1
        count.prob[cols] <- ppois(counts[cols], n.clust.per.col, 
            lower = T)/ppois(n.clust.per.col, n.clust.per.col, 
            lower = T)
        count.prob <- count.prob^count.power["cols"]
    }
    if (iter >= n.iter) {
        cols <- unique(c(cols, names(which(probs >= 0.5))))
    }
    else {
        probs.add <- exp(-(1 - probs)/TEMP)
        probs.drop <- exp(-probs/TEMP)
        col.memb <- attr(ratios, "cnames") %in% cols
        names(col.memb) <- attr(ratios, "cnames")
        sample.probs <- ifelse(col.memb, probs.drop, probs.add)
        if (exists("count.prob")) 
            sample.probs <- sample.probs * count.prob
        n.col.per.clust <- n.clust.per.col/k.clust * length(col.memb)
        balance <- sum(sample.probs[col.memb])/sum(sample.probs[!col.memb]) * 
            n.col.per.clust/sum(col.memb)
        if (is.infinite(balance) || is.na(balance)) 
            balance <- 1
        if (balance >= 1) 
            sample.probs[!col.memb] <- sample.probs[!col.memb] * 
                balance
        else sample.probs[col.memb] <- sample.probs[col.memb]/balance
        allowed.moves <- sample.probs[sample.probs > runif(length(sample.probs))]
        if (length(allowed.moves) > max.change["cols"]) 
            allowed.moves <- sample(allowed.moves, max.change["cols"], 
                prob = sample.probs[names(allowed.moves)])
        cols <- c(cols, names(col.memb[names(allowed.moves)][col.memb[names(allowed.moves)] == 
            FALSE]))
        cols <- cols[!cols %in% names(col.memb[names(allowed.moves)][col.memb[names(allowed.moves)] == 
            TRUE])]
    }
    list(r = rows, c = cols)
}
get.clusterStack <-
function (ks = 1:k.clust, force = F, ...) 
{
    if (!force && !is.null(attr(clusterStack, "iter")) && attr(clusterStack, 
        "iter") == iter) 
        return(clusterStack)
    mc <- get.parallel(length(ks))
    clusterStack <- mc$apply(ks, get.clust, ...)
    tmp <- list.reference(clusterStack, file = sprintf("%s/clusterStack", 
        cmonkey.filename), type = "RDS")
    clusterStack <- tmp
    attr(clusterStack, "iter") <- iter
    clusterStack
}
get.col.scores <-
function (k, for.cols = "all", ratios = ratios[[1]], method = c("new", 
    "orig", "ent")[1], norm.method = c("mean", "all.colVars", 
    "none")[1], ...) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (for.cols[1] == "all") 
        for.cols <- colnames(ratios)
    rows <- rows[rows %in% rownames(ratios)]
    if (length(rows) <= 1) 
        return(rep(NA, length(for.cols)))
    rats <- ratios[rows, for.cols, drop = F]
    row.weights <- if (exists("get.row.weights")) 
        get.row.weights(rows, cols, ratios)
    else NA
    if (method == "orig") {
        if (is.na(row.weights[1])) {
            rats.mn <- colMeans(rats, na.rm = T)
        }
        else {
            rats.mn <- apply(rats, 2, weighted.mean, w = row.weights[rows], 
                na.rm = T)
        }
        rats[, ] <- t(t(rats) - rats.mn)^2
        rats <- colMeans(rats, na.rm = T)
    }
    else if (method == "new") {
        mn <- mean(get.row.scores(k, ratios = ratios, method = row.score.func, 
            ...), na.rm = T)
        rows <- get.rows(k)
        cols <- get.cols(k)
        rats <- -sapply(for.cols, function(cc) {
            if (is.na(row.weights[1])) {
                if (cc %in% cols) 
                  mn/mean(get.row.scores(rows, cols = cols[cols != 
                    cc], ratios = ratios, method = row.score.func, 
                    ...), na.rm = T)
                else mean(get.row.scores(k, cols = c(cols, cc), 
                  ratios = ratios, method = row.score.func, ...), 
                  na.rm = T)/mn
            }
            else {
                if (cc %in% cols) 
                  mn/weighted.mean(get.row.scores(rows, cols = cols[cols != 
                    cc], ratios = ratios, method = row.score.func, 
                    ...), w = row.weights, na.rm = T)
                else weighted.mean(get.row.scores(rows, cols = c(cols, 
                  cc), ratios = ratios, method = row.score.func, 
                  ...), w = row.weights, na.rm = T)/mn
            }
        })
        return(rats)
    }
    var.norm <- 0.99
    if (norm.method == "all.colVars") {
        all.colVars <- attr(ratios, "all.colVars")
        if (!is.null(all.colVars)) 
            var.norm <- all.colVars[for.cols]
    }
    else if (norm.method == "mean") {
        if (!exists("rats.mn")) {
            row.weights <- if (exists("get.row.weights")) 
                get.row.weights(rows, cols, ratios)
            else NA
            rats.tmp <- ratios[rows, for.cols, drop = F]
            if (is.na(row.weights[1])) {
                rats.mn <- colMeans(rats.tmp, na.rm = T)
            }
            else {
                rats.mn <- apply(rats.tmp, 2, weighted.mean, 
                  w = row.weights[rows], na.rm = T)
            }
        }
        var.norm <- abs(rats.mn)
    }
    rats <- rats/(var.norm + 0.01)
    rats
}
get.col.weights <-
function (rows, cols, ratios) 
{
    NA
}
get.cols <-
function (k) 
{
    out <- clusterStack[[k]]$cols
    if (is.null(out) || is.na(out) || length(out) <= 0) 
        out <- attr(ratios, "cnames")
    out
}
get.combined.scores <-
function (quantile.normalize = F) 
{
    r.scores <- row.scores[, ]
    r.scores <- matrix.reference(r.scores)
    if (!quantile.normalize) {
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        rownames(row.memb) <- attr(ratios, "rnames")
        tmp <- r.scores[, ] < -20
        r.scores[, ][tmp] <- min(r.scores[, ][!tmp], na.rm = T)
        rsm <- r.scores[, ][row.memb]
        tmp <- mad(rsm, na.rm = T)
        if (tmp != 0) 
            r.scores[, ] <- (r.scores[, ] - median(rsm, na.rm = T))/tmp
        else {
            tmp <- sd(rsm, na.rm = T)
            if (tmp != 0) 
                r.scores[, ] <- (r.scores[, ] - median(rsm, na.rm = T))/tmp
        }
        rm(tmp, rsm)
    }
    tmp <- r.scores[, ] < -20
    r.scores[, ][tmp] <- min(r.scores[, ][!tmp], na.rm = T)
    rm(tmp)
    r.scores[, ][is.infinite(r.scores[, ])] <- NA
    r.scores[, ][is.na(r.scores[, ])] <- max(r.scores[, ], na.rm = T)
    if (!quantile.normalize && !is.null(mot.scores) || !is.null(net.scores)) 
        rs.quant <- quantile(r.scores[, ], 0.01, na.rm = T)
    if (!is.null(mot.scores)) {
        m.scores <- mot.scores[, ]
    }
    else m.scores <- NULL
    if (!is.null(mot.scores) && !is.null(m.scores)) {
        tmp <- m.scores[, ] < -20
        m.scores[, ][tmp] <- min(m.scores[, ][!tmp], na.rm = T)
        rm(tmp)
        if (!quantile.normalize) {
            m.scores[, ] <- m.scores[, ] - quantile(m.scores[, 
                ], 0.99, na.rm = T)
            m.scores[, ] <- m.scores[, ]/abs(quantile(m.scores[, 
                ], 0.01, na.rm = T)) * abs(rs.quant)
        }
    }
    if (!is.null(net.scores)) {
        n.scores <- net.scores[, ]
        n.scores <- matrix.reference(n.scores)
    }
    else n.scores <- NULL
    if (!is.null(net.scores) && !is.null(n.scores)) {
        n.scores[, ] <- n.scores[, ] - quantile(n.scores[, ], 
            0.99, na.rm = T)
        if (!quantile.normalize) {
            qqq <- abs(quantile(n.scores[, ], 0.01, na.rm = T))
            if (qqq == 0) 
                qqq <- sort(n.scores[, ])[10]
            if (qqq == 0) 
                qqq <- min(n.scores[, ], na.rm = T)
            if (qqq != 0) 
                n.scores[, ] <- n.scores[, ]/qqq * abs(rs.quant)
            rm(qqq)
        }
    }
    if (!is.null(col.scores)) {
        c.scores <- col.scores[, ] * 0
        c.scores <- matrix.reference(c.scores)
        tmp <- c.scores[, ] < -20
        c.scores[, ][tmp] <- min(c.scores[, ][!tmp], na.rm = T)
        rm(tmp)
    }
    else c.scores <- NULL
    new.weights <- c(row = row.scaling[iter], mot = mot.scaling[iter], 
        net = net.scaling[iter])
    if (pareto.adjust.scalings && iter > 51) {
        new.weights <- pareto.adjust.weights()
        if (iter %in% stats.iters) 
            cat("New weights:", new.weights, "\n")
    }
    if (quantile.normalize) {
        tmp <- list(row = r.scores, mot = m.scores, net = n.scores)
        if (sum(sapply(tmp, function(i) !is.null(i))) > 1) {
            wts <- new.weights[!sapply(tmp, is.null)]
            tmp <- quantile.normalize.scores(tmp, weights = wts)
            if (!is.null(r.scores)) 
                r.scores[, ] <- tmp$row[, ]
            if (!is.null(m.scores)) 
                m.scores[, ] <- tmp$mot[, ]
            if (!is.null(n.scores)) 
                n.scores[, ] <- tmp$net[, ]
        }
        rm(tmp)
    }
    if (new.weights["row"] != 1) 
        r.scores[, ] <- r.scores[, ] * new.weights["row"]
    if (!is.null(m.scores)) {
        tmp <- !is.na(m.scores[, ])
        r.scores[, ][tmp] <- r.scores[, ][tmp] + m.scores[, ][tmp] * 
            new.weights["mot"]
    }
    if (!is.null(n.scores)) {
        tmp <- !is.na(n.scores[, ])
        r.scores[, ][tmp] <- r.scores[, ][tmp] + n.scores[, ][tmp] * 
            new.weights["net"]
    }
    r.scores <- matrix.reference(r.scores)
    c.scores <- matrix.reference(c.scores)
    invisible(list(r = r.scores, c = c.scores, scalings = new.weights))
}
get.condition.groups <-
function (conds = attr(ratios, "cnames")) 
{
    tmp <- sapply(conds, function(i) strsplit(i, "__")[[1]][1])
    tmp2 <- as.integer(as.factor(tmp))
    names(tmp2) <- names(tmp)
    tmp2
}
get.conditions <-
function (biclusters) 
{
    if (!missing(biclusters)) {
        inds <- as.integer(sapply(strsplit(biclusters, "_"), 
            "[", 2))
        out <- lapply(inds, function(i) e$clusterStack[[i]]$cols)
        names(out) <- biclusters
    }
    out
}
get.density.scores <-
function (ks = 1:k.clust, r.scores, c.scores, plot = "none", 
    bw.scale = function(nr) exp(-nr/10) * 10) 
{
    rr <- attr(ratios, "rnames")
    rs <- r.scores
    bw.r <- max(diff(range(rs[, ], na.rm = T))/100, 0.001)
    get.rr.scores <- function(k) {
        rows <- get.rows(k)
        cols <- get.cols(k)
        rsk <- rs[, k, drop = T]
        if (length(rows) > 0 && length(cols) > 0) {
            bw <- bw.r * bw.scale(length(rows))
            d <- density(rsk[rows], na.rm = T, bw = bw, adjust = 2, 
                from = min(rsk, na.rm = T) - 1, to = max(rsk, 
                  na.rm = T) + 1, n = 256)
            p <- approx(d$x, rev(cumsum(rev(d$y))), rsk)$y
            if ("rows" %in% plot) {
                h = hist(rsk, breaks = 50, main = NULL, xlab = "Combined scores")
                tmp.scale <- round(attr(ratios, "nrow")/length(rows)/4)
                hist(rep(rsk[rows], tmp.scale), breaks = h$breaks, 
                  col = "red", border = "red", add = T)
                hist(rsk, breaks = h$breaks, add = T)
                lines(d$x, d$y/max(d$y, na.rm = T) * attr(ratios, 
                  "nrow")/50, col = "blue")
                lines(sort(rsk), p[order(rsk)]/max(p, na.rm = T) * 
                  attr(ratios, "nrow")/50, col = "green")
            }
        }
        else {
            p <- rep(1, attr(ratios, "nrow"))
        }
        return(p/sum(p, na.rm = T))
    }
    rr.scores <- row.scores[, ] * 0
    rr.scores <- matrix.reference(rr.scores)
    mc <- get.parallel(length(ks))
    rr.scores[, ] <- do.call(cbind, mc$apply(ks, get.rr.scores))
    rr.scores[, ][is.infinite(rr.scores[, ])] <- NA
    cc.scores <- NULL
    if (!is.null(col.scores)) {
        cs <- c.scores
        bw.c <- max(diff(range(cs[, ], na.rm = T))/100, 0.001)
        get.cc.scores <- function(k) {
            cols <- get.cols(k)
            rows <- get.rows(k)
            csk <- cs[, k, drop = T]
            if (length(cols) > 0 && length(rows) > 0 && !all(is.na(csk[cols])) && 
                !all(is.infinite(csk[cols])) & !all(csk[cols][!is.na(csk[cols])] == 
                csk[cols[!is.na(csk[cols])][1]])) {
                d <- density(csk[cols], na.rm = T, from = min(csk, 
                  na.rm = T) - 1, to = max(csk, na.rm = T) + 
                  1, bw = bw.c, adjust = 2, n = 256)
                p <- approx(d$x, rev(cumsum(rev(d$y))), csk)$y
                if ("cols" %in% plot) {
                  h = hist(csk, breaks = 50, main = NULL, xlab = "Combined scores")
                  tmp.scale <- round(attr(ratios, "ncol")/length(cols)/4) + 
                    1
                  hist(rep(csk[cols], tmp.scale), breaks = h$breaks, 
                    col = "red", border = "red", add = T)
                  hist(csk, breaks = h$breaks, add = T)
                  lines(d$x, d$y/max(d$y, na.rm = T) * attr(ratios, 
                    "ncol")/50, col = "blue")
                  lines(sort(csk), p[order(csk)]/max(p, na.rm = T) * 
                    attr(ratios, "ncol")/50, col = "green")
                }
            }
            else {
                p <- rep(1, attr(ratios, "ncol"))
            }
            return(p/sum(p, na.rm = T))
        }
        cc.scores <- col.scores[, ] * 0
        cc.scores <- matrix.reference(cc.scores)
        if (!is.null(c.scores)) {
            cc.scores[, ] <- do.call(cbind, mc$apply(ks, get.cc.scores))
            cc.scores[, ][is.infinite(cc.scores)] <- NA
        }
    }
    invisible(list(r = rr.scores, c = cc.scores))
}
get.dup.seqs <-
function (seqs) 
{
    out <- duplicated(seqs)
    names(out) <- names(seqs)
    out
}
get.expression <-
function (genes, biclusters) 
{
    if (!missing(genes)) {
        out <- lapply(genes, function(g) e$ratios$ratios[g, ])
        names(out) <- genes
    }
    if (!missing(biclusters)) {
        out <- lapply(biclusters, function(bic) {
            print(bic)
            e$ratios$ratios[get.genes(biclusters = bic)[[1]], 
                get.conditions(biclusters = bic)[[1]]]
        })
        names(out) <- biclusters
    }
    out
}
get.expressions <-
function (genes, biclusters) 
{
    if (!missing(genes)) {
        out <- lapply(genes, function(g) e$ratios$ratios[g, ])
        names(out) <- genes
    }
    if (!missing(biclusters)) {
        out <- lapply(biclusters, function(bic) {
            print(bic)
            e$ratios$ratios[get.genes(biclusters = bic)[[1]], 
                get.conditions(biclusters = bic)[[1]]]
        })
        names(out) <- biclusters
    }
    out
}
get.gene.coords <-
function (rows, op.shift = T, op.table = genome.info$operons, 
    ...) 
{
    rows <- unique(rows)
    syns <- get.synonyms(rows, ...)
    tab <- genome.info$feature.tab
    ids <- lapply(syns, function(s) s[s %in% tab$id])
    if (all(sapply(ids, length) < 1)) {
        warning("Could not find gene start/stop for any input genes", 
            call. = F)
        return(NULL)
    }
    if (any(sapply(ids, length) < 1)) 
        warning("Could not find gene start/stop for all input genes", 
            call. = F)
    ids <- ids[sapply(ids, length) >= 1]
    ids <- sapply(ids, "[", 1)
    ids <- data.frame(id = ids, names = names(ids))
    coos <- NULL
    if (op.shift) {
        if (attr(op.table, "source") == "rsat") {
            ops <- merge(ids, op.table, by.x = "id", by.y = "query", 
                all = F)
            ops2 <- ops[order(ops$lead), ]
            coos <- merge(ops, tab, by.x = "lead", by.y = "name", 
                all = F)[, c("id.x", "names", "contig", "strand", 
                "start_pos", "end_pos")]
        }
        else if (attr(op.table, "source") == "microbes.online") {
            ops <- NULL
            if (!any(ids$names %in% op.table$gene)) {
                ids2 <- lapply(syns, function(s) s[s %in% op.table$gene])
                if (all(sapply(ids2, length) < 1)) {
                  warning("Could not find operon info for any input genes", 
                    call. = F)
                }
                else {
                  if (any(sapply(ids2, length) < 1)) 
                    warning("Could not find operon info for all input genes", 
                      call. = F)
                  ids2 <- ids2[sapply(ids2, length) >= 1]
                  ids2 <- sapply(ids2, "[", 1)
                  ids2 <- data.frame(id = ids2, names = names(ids2))
                  ops <- merge(ids2, op.table, by.x = "id", by.y = "gene", 
                    all.x = T)
                }
            }
            if (is.null(ops)) 
                ops <- merge(ids, op.table, by.x = "names", by.y = "gene", 
                  all.x = T)
            if (any(is.na(ops$head))) {
                head <- as.character(ops$head)
                head[is.na(head)] <- as.character(ops$names[is.na(head)])
                ops$head <- as.factor(head)
            }
            head.syns <- get.synonyms(unique(as.character(ops$head)))
            head.ids <- lapply(head.syns, function(s) s[s %in% 
                tab$id])
            head.ids <- head.ids[sapply(head.ids, length) >= 
                1]
            head.ids <- data.frame(id = sapply(head.ids, "[", 
                1), names = names(head.ids))
            ops2 <- merge(ops, head.ids, by.x = "head", by.y = "names", 
                all.x = T)
            coos <- merge(ops2, tab, by.x = "id.y", by.y = "id", 
                all.x = T)[, c("id.x", "names", "contig", "strand", 
                "start_pos", "end_pos")]
        }
    }
    else {
        coos <- merge(ids, tab, by = "id")[, c("id", "names", 
            "contig", "strand", "start_pos", "end_pos")]
    }
    colnames(coos)[1] <- "id"
    if (is.factor(coos$start_pos)) 
        coos$start_pos <- as.numeric(levels(coos$start_pos))[coos$start_pos]
    if (is.factor(coos$end_pos)) 
        coos$end_pos <- as.numeric(levels(coos$end_pos))[coos$end_pos]
    coos[!duplicated(coos[, 1:4]), ]
}
get.genes <-
function (biclusters, motifs) 
{
    if (!missing(biclusters)) {
        inds <- as.integer(sapply(strsplit(biclusters, "_"), 
            "[", 2))
        out <- lapply(inds, function(i) e$clusterStack[[i]]$rows)
        names(out) <- biclusters
    }
    else if (!missing(motifs)) {
        out <- lapply(motifs, function(m) {
            tmp <- get.motif.info(motif = m)[[1]]
            if (is.null(tmp)) 
                return(NULL)
            if ("pvals" %in% colnames(tmp$mast)) 
                unique(as.character(subset(tmp$mast, pvals <= 
                  0.05 & abs(mots) == as.integer(strsplit(m, 
                  "_")[[1]][3]), gene, drop = T)))
            else unique(as.character(subset(tmp$mast, p.value <= 
                0.05, gene, drop = T)))
        })
        names(out) <- motifs
    }
    out
}
get.genome.info <-
function (fetch.upstream = F, fetch.predicted.operons = "rsat") 
{
    rsat.url <- rsat.urls[1]
    feature.tab <- feature.names <- genome.seqs <- operons <- org.id <- synonyms <- NULL
    genome.loc <- paste(rsat.url, "/data/genomes/", rsat.species, 
        "/genome/", sep = "")
    fname <- paste("data/", rsat.species, "/organism_names.tab", 
        sep = "")
    err <- dlf(fname, paste(genome.loc, "/organism_names.tab", 
        sep = ""))
    if (class(err) == "try-error") {
        tmp.url <- paste(rsat.url, "/data/genomes/", rsat.species, 
            "_EnsEMBL/genome/organism_names.tab", sep = "")
        err <- dlf(fname, tmp.url)
        if (class(err) != "try-error") 
            genome.loc <- paste(rsat.url, "/data/genomes/", rsat.species, 
                "_EnsEMBL/genome/", sep = "")
    }
    if (!file.exists(fname) || file.info(fname)$size <= 0) 
        stop(paste("Genome info for", rsat.species, "does not exist. Please check", 
            genome.loc, "and let me know if I am wrong"))
    nskip <- sum(substr(readLines(gzfile(fname), n = 20), 1, 
        2) == "--" | readLines(gzfile(fname), n = 20) == "")
    org.id <- read.delim(gzfile(fname), head = F, as.is = T, 
        skip = nskip)
    if (!exists("taxon.id") || is.na(taxon.id) || is.null(taxon.id)) 
        taxon.id <- org.id$V1[1]
    cat("Organism taxon id:", taxon.id, "\n")
    closeAllConnections()
    if (!no.genome.info && !all(grepl("file=", names(mot.weights)))) {
        fname <- paste("data/", rsat.species, "/feature.tab", 
            sep = "")
        use.cds <- FALSE
        err <- dlf(fname, paste(genome.loc, "feature.tab", sep = ""), 
            paste("Fetching genome annotation data from RSAT", 
                rsat.url, "..."))
        if (class(err) == "try-error") {
            err <- dlf(fname, paste(genome.loc, "cds.tab", sep = ""))
            use.cds <- TRUE
        }
        cat("Loading genome annotation data...\n")
        head <- readLines(gzfile(fname), n = 30)
        nskip <- length(grep("^--", head))
        feature.tab <- read.delim(gzfile(fname), skip = nskip, 
            head = F, comment = "", as.is = F)
        closeAllConnections()
        head <- strsplit(gsub("^-- ", "", head[grep("^-- id", 
            head, perl = T)], perl = T), "\t")[[1]]
        colnames(feature.tab) <- head[1:ncol(feature.tab)]
        fname <- paste("data/", rsat.species, "/feature_names.tab", 
            sep = "")
        err <- dlf(fname, paste(genome.loc, if (!use.cds) 
            "feature_names.tab"
        else "cds_names.tab", sep = ""))
        nskip <- sum(substr(readLines(gzfile(fname), n = 20), 
            1, 2) == "--")
        closeAllConnections()
        feature.names <- read.delim(gzfile(fname), head = F, 
            as.is = T, skip = nskip, row.names = NULL, comment = "")
        closeAllConnections()
        colnames(feature.names) <- c("id", "names", "type")
        feature.names <- unique(feature.names)
        chroms <- unique(as.character(feature.tab$contig))
        chroms <- chroms[!is.na(chroms) & chroms != ""]
        if (!is.na(mot.iters[1])) {
            genome.seqs <- list()
            for (i in chroms) {
                cat("Loading genome sequence, chromosome", i, 
                  "\n")
                fname <- paste("data/", rsat.species, "/", i, 
                  ".raw", sep = "")
                err <- dlf(fname, paste(genome.loc, i, ".raw", 
                  sep = ""))
                if (class(err) == "try-error") {
                  ii <- gsub(":", "_", i, fixed = T)
                  err <- dlf(fname, paste(genome.loc, ii, ".raw", 
                    sep = ""))
                  if (class(err) == "try-error") {
                    err <- dlf(fname, paste(genome.loc, gsub(".[0-9]$", 
                      "", i), ".raw", sep = ""))
                    if (class(err) == "try-error") 
                      cat("ERROR reading genome sequence", i, 
                        "\n")
                    else fname <- paste("data/", rsat.species, 
                      "/", gsub(".[0-9]$", "", i), ".raw", sep = "")
                  }
                  else fname <- paste("data/", rsat.species, 
                    "/", ii, ".raw", sep = "")
                }
                out <- try(readLines(gzfile(fname)), silent = T)
                if (class(out) == "try-error" || length(out) == 
                  0 || is.na(out) || out == "" || out == "CHARACTER(0)") 
                  out <- try(readLines(fname), silent = T)
                if (class(out) == "try-error" || length(out) == 
                  0 || is.na(out) || out == "" || out == "CHARACTER(0)") {
                  cat("ERROR reading genome sequence", i, "\n")
                  next
                }
                out <- toupper(out)
                if (FALSE && big.memory == TRUE) {
                  rf <- as.factor(strsplit(out, character(0))[[1]])
                  require(ff)
                  frf <- ff(rf, vmode = "quad", levels = levels(rf), 
                    filename = sprintf("./%s/%s.genome.ff", cmonkey.filename, 
                      i))
                  out <- frf
                }
                genome.seqs[[i]] <- out
            }
            if (length(genome.seqs) != length(chroms)) {
                cat("WARNING: Could not read sequence for chromosomes", 
                  chroms[!chroms %in% names(genome.seqs)], "\n")
                feature.tab <- subset(feature.tab, contig %in% 
                  names(genome.seqs))
            }
            if (length(genome.seqs) <= 0) 
                genome.seqs <- NULL
        }
        if (!is.na(mot.iters[1]) && fetch.upstream) {
            fname <- paste("data/", rsat.species, "/upstream-noorf.fasta.gz", 
                sep = "")
            err <- dlf(fname, paste(genome.loc, rsat.species, 
                "_upstream-noorf.fasta.gz", sep = ""), "Fetching upstream sequences from RSAT...")
            upstream.noorf <- readLines(gzfile(fname))
            fname <- paste("data/", rsat.species, "/upstream.fasta.gz", 
                sep = "")
            err <- dlf(fname, paste(genome.loc, rsat.species, 
                "_upstream.fasta.gz", sep = ""))
            upstream <- readLines(gzfile(fname))
        }
    }
    synonyms <- NULL
    if (exists("synonym.thesaurus")) {
        cat("Loading synonym information from custom thesaurus...\n")
        tmp <- read.csv(gzfile(synonym.thesaurus), header = F)
        synonyms <- strsplit(toupper(as.character(tmp[, 2])), 
            ";")
        names(synonyms) <- toupper(as.character(tmp[, 1]))
        rm(tmp)
    }
    else if (exists("ratios") && !is.null(ratios)) {
        cat("Gathering all \"standard\" orf names and other synonyms for all probe names...\n")
        synonyms <- get.synonyms(attr(ratios, "rnames"), feature.names, 
            verbose = T)
    }
    if (!is.null(synonyms)) {
        cat("Mean number of synonyms per probe:", mean(sapply(synonyms, 
            length), na.rm = T), "\n")
        is.bad <- sapply(names(synonyms), function(i) length(synonyms[[i]]) == 
            0 || substr(synonyms[[i]][1], 1, 5) == "Error")
        if (sum(is.bad) > 0) {
            cat("These", sum(is.bad), "probe names have no matching ORF annotation:\n")
            print(names(which(is.bad)))
        }
        rm(is.bad)
    }
    closeAllConnections()
    invisible(list(species = rsat.species, genome.seqs = genome.seqs, 
        feature.tab = feature.tab, feature.names = feature.names, 
        org.id = org.id, taxon.id = taxon.id, synonyms = synonyms))
}
get.long.names <-
function (k, shorter = F) 
{
    if (is.numeric(k[1])) {
        rows <- get.rows(k)
    }
    else {
        rows <- k
    }
    if (is.null(genome.info$feature.tab)) {
        out <- rep("", length(rows))
        names(out) <- rows
        return(rows)
    }
    ids <- get.synonyms(rows)
    mc <- list(apply = lapply)
    if (!shorter) 
        desc <- mc$apply(ids, function(i) subset(genome.info$feature.tab, 
            id %in% i, select = c("id", "description")))
    else {
        desc <- mc$apply(ids, function(i) subset(genome.info$feature.tab, 
            id %in% i, select = c("id", "name", "description")))
        for (i in 1:length(desc)) if (length(desc[[i]]$name) > 
            0 && desc[[i]]$name %in% rows) {
            if (grepl("(", desc[[i]]$description, fixed = T)) 
                desc[[i]]$name <- strsplit(as.character(desc[[i]]$description), 
                  "[()]", perl = T)[[1]][2]
        }
    }
    out <- sapply(desc, function(i) as.character(i[1, 2]))
    out <- out[rows]
    names(out) <- rows
    out[is.na(out) | out == names(out)] <- ""
    out
}
get.mast.pvals <-
function (mast.output, in.genes = NULL) 
{
    space.pad <- function(lines, length) {
        nc <- nchar(lines)
        nc[nc >= length] <- 0
        spaces <- sapply(1:length(lines), function(i) paste(rep(" ", 
            length - nc[i]), sep = "", collapse = ""))
        paste(lines, spaces)
    }
    out <- list()
    start <- grep("SECTION III: ANNOTATED SEQUENCES", mast.output)
    if (length(start) == 0 || is.na(start)) 
        return(out)
    end <- grep("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", 
        mast.output[(start + 3):length(mast.output)]) + start + 
        3
    line.starts <- grep("LENGTH = ", mast.output[(start + 2):(start + 
        1 + end)]) + start + 1
    if (is.null(line.starts) || length(line.starts) == 0) 
        return(out)
    for (i in 1:length(line.starts)) {
        l <- line.starts[i]
        gene <- mast.output[l - 2]
        if (is.null(gene) || is.na(gene) || (!is.null(in.genes) && 
            !(gene %in% in.genes))) 
            next
        l.next <- line.starts[i + 1] - 2
        if (i >= length(line.starts)) 
            l.next <- end
        if (l.next - l <= 5) 
            next
        submast <- mast.output[l:(l.next - 1)]
        l.start <- which(submast == "")[1] + 1
        if (submast[l.start] == "") 
            l.start <- l.start + 1
        q <- list()
        for (i in 1:6) q[[i]] <- space.pad(submast[seq((l.start + 
            i - 1), length(submast), by = 6)], 80)
        seq.starts <- as.integer(sapply(strsplit(q[[5]], " "), 
            "[", 1))
        char.skip <- which(strsplit(q[[5]][1], "")[[1]] %in% 
            c("G", "A", "T", "C", "N", "X"))[1]
        mots <- unlist(strsplit(gsub("[\\[\\]\\<\\>]", "", paste(substr(q[[1]], 
            char.skip, 80), collapse = ""), perl = T), "\\s+", 
            perl = T))
        mots <- as.integer(mots[!is.na(as.integer(mots))])
        mots <- mots[!is.na(mots)]
        p.vals <- strsplit(paste(substr(q[[2]], char.skip, 80), 
            collapse = ""), "\\s+")[[1]]
        p.vals <- as.numeric(p.vals[!is.na(as.numeric(p.vals))])
        posns <- integer()
        for (i in 1:length(q[[1]])) {
            posns <- c(posns, which(strsplit(substr(q[[1]][i], 
                char.skip, 80), "")[[1]] %in% c("[", "<")) + 
                seq.starts[i])
        }
        out[[gene]] <- list(pvals = p.vals, mots = mots, posns = posns)
    }
    return(out)
}
get.motif.cluster.info <-
function (motif.clusters) 
{
    if (!missing(motif.clusters)) {
        out <- tt.out2[as.integer(gsub("MOTC_", "", motif.clusters))]
        names(out) <- motif.clusters
    }
    out
}
get.motif.clusters <-
function (motifs) 
{
    if (!missing(motifs)) {
        tmp <- get.motifs(motif.clust = paste("MOTC", 1:length(tt.out2), 
            sep = "_"))
        out <- lapply(motifs, function(m) names(which(sapply(tmp, 
            function(i) m %in% i))))
        names(out) <- motifs
    }
    out
}
get.motif.info <-
function (motifs) 
{
    if (!missing(motifs)) {
        inds <- strsplit(motifs, "_")
        out <- lapply(inds, function(i) {
            ms <- e$meme.scores[[1]][[as.integer(i[2])]]
            if (is.null(ms$meme.out) || length(ms$meme.out) < 
                as.integer(i[3])) 
                return(NULL)
            out <- ms$meme.out[[as.integer(i[3])]]
            tmp <- ms$pv.ev[[2]]
            if (!is.null(tmp)) 
                tmp <- subset(tmp, abs(mots) == as.integer(i[3]))
            else tmp <- out$posns
            out$mast <- tmp
            out
        })
        names(out) <- motifs
    }
    out
}
get.motif.positions <-
function (motifs = "ALL", seqs = e$genome.info$all.upstream.seqs[[1]], 
    seq.type = names(e$meme.scores)[1], counts = T, verbose = F) 
{
    if (!exists("pssm.scans")) 
        pssm.scans <- get.pssm.scans(motifs, seqs, seq.type)
    p.scans <- pssm.scans
    p.scans$mots <- abs(p.scans$mots)
    gc()
    p.scans <- p.scans[, c("bic", "mots", "gene", "posns")]
    gc()
    require(data.table)
    p.scans <- as.data.table(p.scans)
    gc()
    setkey(p.scans, bic, mots, gene, posns)
    if (motifs != "ALL") {
        mots <- strsplit(gsub("MOT_", "", motifs), "_")
        bi <- as.integer(sapply(mots, "[", 1))
        mo <- as.integer(sapply(mots, "[", 2))
        scans <- p.scans[J(bi, mo, names(seqs))]
    }
    else {
        scans <- p.scans
    }
    scans <- scans[!is.na(posns)]
    if (!counts) 
        return(invisible(scans))
    ms <- e$meme.scores[[seq.type]]
    counts <- integer()
    if (nrow(scans) > 0) {
        for (i in 1:nrow(scans)) {
            if (verbose && i%%1000 == 0) 
                cat(i, nrow(scans), length(counts), "\n")
            k <- scans$bic[i]
            mot <- scans$mots[i]
            width <- ms[[k]]$meme.out[[mot]]$width
            if (width > 0) {
                posn <- scans$posns[i]
                inds <- (posn - 1):(posn - 2 + width)
                if (length(counts) < inds[length(inds)]) {
                  counts[inds[length(inds)]] <- 0
                  counts[is.na(counts)] <- 0
                }
                counts[inds] <- counts[inds] + 1
            }
        }
    }
    invisible(counts)
}
get.motif.scores <-
function (k, out.ms, for.rows = "all") 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (for.rows[1] == "all") 
        for.rows <- attr(ratios, "rnames")
    if (length(rows) <= 1 || is.null(out.ms$pv.ev)) 
        return(rep(NA, length(for.rows)))
    m.scores <- rep(NA, attr(ratios, "nrow"))
    names(m.scores) <- attr(ratios, "rnames")
    tmp <- out.ms$pv.ev[[1]][, "p.value"]
    names(tmp) <- rownames(out.ms$pv.ev[[1]])
    m.scores[names(tmp)] <- tmp
    m.scores <- log(m.scores)
    m.scores <- m.scores[for.rows]
    return(m.scores)
}
get.motifs <-
function (biclusters, motif.clusters, genes, positions, window = 10) 
{
    if (!missing(biclusters)) {
        out <- lapply(biclusters, function(bic) paste(gsub("BIC_", 
            "MOT_", bic), 1:length(e$clusterStack[[as.integer(gsub("BIC_", 
            "", bic))]]$e.val), sep = "_"))
        names(out) <- biclusters
    }
    if (!missing(motif.clusters)) {
        tmp <- as.integer(sapply(strsplit(motif.clusters, "_"), 
            "[", 2))
        out <- lapply(tmp, function(i) paste("MOT", attr(tt.out2[[i]], 
            "mot.names"), sep = "_"))
        names(out) <- motif.clusters
    }
    if (!missing(genes)) {
        out <- lapply(genes, function(g) {
            unlist(lapply(paste("MOT", c(paste(1:e$k.clust, 1, 
                sep = "_"), paste(1:e$k.clust, 2, sep = "_")), 
                sep = "_"), function(m) {
                tmp <- get.motif.info(motif = m)[[1]]
                if (is.null(tmp)) 
                  return(NULL)
                if ("pvals" %in% colnames(tmp$mast)) 
                  tmp <- unique(as.character(subset(tmp$mast, 
                    pvals <= 0.05 & abs(mots) == as.integer(strsplit(m, 
                      "_")[[1]][3]), gene, drop = T)))
                else tmp <- unique(as.character(subset(tmp$mast, 
                  p.value <= 0.05, gene, drop = T)))
                if (g %in% tmp) 
                  return(m)
                return(NULL)
            }))
        })
        names(out) <- genes
    }
    if (!missing(positions)) {
        if (!exists("pssm.scans")) 
            pssm.scans <- get.motif.positions(NULL)
        if (positions[1] == "choose") {
            pn <- names(positions)
            positions <- locator(1, "p")$x
            names(positions) <- pn
        }
        else if (grepl("[, ]", positions[1], perl = T)) {
            positions <- strsplit(positions, "[, ]", perl = T)[[1]]
        }
        if (length(positions) == 3) {
            tmp <- as.integer(positions[2:3])
            names(tmp)[1] <- positions[1]
            positions <- tmp
        }
        else if (is.character(positions) && length(positions) == 
            2) {
            positions <- as.integer(positions)
        }
        if (is.null(names(positions)[1])) 
            names(positions)[1] <- "Chr"
        if (length(positions) == 1) {
            tmp <- c(positions[1] - window, positions[1] + window)
            names(tmp)[1] <- names(positions)[1]
            positions <- tmp
        }
        print(positions)
        chr <- names(positions)[1]
        scans <- subset(pssm.scans, gene == chr & posns %betw% 
            positions)
        out <- paste("MOT", scans$bic, scans$mots, sep = "_")
    }
    out
}
get.network.scores <-
function (k, net = networks$string, for.rows = "all", p1.col = "protein1", 
    p2.col = "protein2", score.col = "combined_score", combine.func = sum) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (for.rows[1] == "all") 
        for.rows <- attr(ratios, "rnames")
    if (length(rows) < 1) 
        return(rep(NA, length(for.rows)))
    cons <- net[as.character(net[[p1.col]]) %in% rows, c(p2.col, 
        score.col), drop = F]
    if (is.null(cons) || nrow(cons) <= 0) 
        return(rep(NA, length(for.rows)))
    cons <- cons[as.character(cons[[p2.col]]) %in% for.rows, 
        , drop = F]
    if (is.null(cons) || nrow(cons) <= 0) 
        return(rep(NA, length(for.rows)))
    tmp <- tapply(as.numeric(cons[[score.col]]), as.character(cons[[p2.col]]), 
        combine.func, na.rm = T)/length(rows)
    scores <- rep(NA, length(for.rows))
    names(scores) <- for.rows
    scores[names(tmp)] <- tmp
    return(-log(scores + 1))
}
get.old.scores.matrices <-
function (ks = 1:k.clust) 
{
    mc <- get.parallel(length(ks))
    all.scores <- mc$apply(ks, get.all.scores, return.scores = T, 
        densities = F, verbose = F, members = F, force.motif = F, 
        allow.motif = F)
    row.scores <- sapply(all.scores, function(i) i$r[, "ratios.ratios"])
    mot.scores <- lapply(all.scores, function(i) i$r[, grep("^motif\\.", 
        colnames(i$r)), drop = F])
    if (!is.null(unlist(mot.scores))) {
        mot.scores <- lapply(colnames(mot.scores[[1]]), function(i) sapply(mot.scores, 
            function(j) j[, i]))
        if (length(mot.scores) > 1) 
            for (m in 2:length(mot.scores)) mot.scores[[1]] <- mot.scores[[1]] + 
                mot.scores[[m]]
        mot.scores <- mot.scores[[1]]
    }
    else rm(mot.scores)
    net.scores <- lapply(all.scores, function(i) i$r[, grep("^network\\.", 
        colnames(i$r)), drop = F])
    if (!is.null(unlist(net.scores))) {
        net.scores <- lapply(colnames(net.scores[[1]]), function(i) sapply(net.scores, 
            function(j) j[, i]))
        if (length(net.scores) > 1) 
            for (m in 2:length(net.scores)) net.scores[[1]] <- net.scores[[1]] + 
                net.scores[[m]]
        net.scores <- net.scores[[1]]
    }
    else rm(net.scores)
    col.scores <- sapply(all.scores, function(i) i$c[, "ratios.ratios"])
    list(r = row.scores, m = mot.scores, n = net.scores, c = col.scores)
}
get.operon.predictions <-
function (fetch.predicted.operons = "microbes.online", org.id = genome.info$org.id$V1[1]) 
{
    operons <- NULL
    if (fetch.predicted.operons == "rsat") {
        rsat.url <- rsat.urls[1]
        cat("Using operon predictions from RSAT...\n")
        fname <- paste("data/", rsat.species, "/rsat_operon_predictions.html", 
            sep = "")
        err <- dlf(fname, paste(rsat.url, "/infer-operons.cgi?organism=", 
            rsat.species, "&genes=all&return_leader=on&return_operon=on&return_query=on&", 
            "output=display&dist_thr=55", sep = ""))
        operons <- readLines(gzfile(fname))
        start <- which(operons == "<INPUT type=\"hidden\" NAME=\"gene_selection\" VALUE=\"#lead\toperon\tquery") + 
            1
        end <- which(operons == "<INPUT type=\"hidden\" NAME=\"feattype\" VALUE=\"\">") - 
            2
        operons <- do.call(rbind, strsplit(operons[start:end], 
            "\t+", perl = T))
        colnames(operons) <- c("lead", "operon", "query")
        operons <- as.data.frame(operons)
    }
    else if (fetch.predicted.operons == "microbes.online") {
        cat("Using operon predictions from MicrobesOnline...\n")
        fname <- paste("data/", rsat.species, "/microbesonline_operons_gnc", 
            org.id, ".named", sep = "")
        err <- dlf(fname, paste("http://www.microbesonline.org/operons/gnc", 
            org.id, ".named", sep = ""))
        if (org.id != taxon.id && (!file.exists(fname) || file.info(fname)$size == 
            0)) {
            fname <- paste("data/", rsat.species, "/microbesonline_operons_gnc", 
                taxon.id, ".named", sep = "")
            err <- dlf(fname, paste("http://www.microbesonline.org/operons/gnc", 
                taxon.id, ".named", sep = ""))
        }
        if (file.exists(fname)) 
            cat("Succesfully fetched operon predictions. Parsing...\n")
        ops <- read.delim(gzfile(fname))
        ops2 <- subset(ops, bOp == "TRUE" & SysName1 != "" & 
            SysName2 != "")
        gns <- sort(unique(c(as.character(ops2$SysName1), as.character(ops2$SysName2))))
        gns <- gns[gns != ""]
        sn1 <- as.character(ops2$SysName1)
        sn1[sn1 == "" | is.na(sn1)] <- as.character(ops2$Name1)[sn1 == 
            "" | is.na(sn1)]
        sn2 <- as.character(ops2$SysName2)
        sn2[sn2 == "" | is.na(sn2)] <- as.character(ops2$Name2)[sn2 == 
            "" | is.na(sn2)]
        operons <- list(0)
        for (i in 1:length(sn1)) {
            sn1i <- sn1[i]
            found <- which(sapply(operons, function(j) sn1i %in% 
                j))
            if (length(found) > 0) 
                operons[[found[1]]] <- c(operons[[found[1]]], 
                  sn2[i])
            else operons[[length(operons) + 1]] <- c(sn1i, sn2[i])
        }
        operons <- operons[-1]
        search.names <- c(gns, as.character(genome.info$feature.names$id))
        if (exists("ratios")) 
            search.names <- c(attr(ratios, "rnames"), search.names)
        mc <- get.parallel(length(operons))
        nms <- mc$apply(1:length(operons), function(i) {
            s <- get.synonyms(operons[[i]])
            s <- lapply(s, function(i) i[i %in% search.names])
            ids <- unlist(lapply(s, function(i) i[i %in% genome.info$feature.names$id][1]))
            if (length(ids) <= 0) {
                warning(paste("No genome annotation for any genes in operon #", 
                  i, " -- don't know what to do!", call. = F))
                return("")
            }
            ids[is.na(ids)] <- names(ids)[is.na(ids)]
            vngs <- unlist(lapply(s, function(i) {
                out <- i[!i %in% genome.info$feature.names$id]
                if (length(out) <= 0) 
                  out <- i[i %in% search.names]
                if (length(out) <= 0) 
                  out <- i[genome.info$feature.names$id == i & 
                    genome.info$feature.names$id == "primary"]
                if (length(out) <= 0) 
                  out <- i
                out
            }))
            coos <- get.gene.coords(ids, op.shift = F)
            vngs <- vngs[ids %in% coos$names]
            if (is.null(coos) || nrow(coos) <= 0) {
                warning(paste("No genome annotation for any genes in operon #", 
                  i, " -- don't know what to do!", call. = F))
                return("")
            }
            if (mean(as.character(coos$strand) == "D") > 0.6) 
                head <- vngs[which.min(coos$start_pos)]
            else if (mean(as.character(coos$strand) == "R") > 
                0.6) 
                head <- vngs[which.max(coos$end_pos)]
            else {
                head <- ""
                warning(paste("About 50% of operon #", i, "are on opposite strands -- don't know what to do!", 
                  call. = F))
            }
            head
        })
        names(operons) <- unlist(nms)
        operons <- operons[names(operons) != ""]
        operons <- do.call(rbind, lapply(names(operons), function(h) data.frame(head = h, 
            gene = operons[[h]])))
        operons <- subset(operons, head != "")
    }
    if (!is.null(operons)) 
        attr(operons, "source") <- fetch.predicted.operons
    closeAllConnections()
    operons
}
get.parallel <-
function (X = k.clust, verbose = F, para.cores = get("parallel.cores")) 
{
    if ((exists("DEBUG") && !is.function(DEBUG) && DEBUG == TRUE) || 
        is.na(para.cores) || (is.logical(para.cores) && para.cores == 
        FALSE) || (is.numeric(para.cores) && para.cores <= 1)) {
        out <- list(mc = FALSE, par = para.cores, apply = lapply)
        if (verbose) 
            cat("NOT PARALLELIZING\n")
    }
    else {
        try(has.multi <- require(multicore, quietly = T), silent = T)
        if (!has.multi || (has.multi && multicore:::isChild())) {
            out <- list(mc = FALSE, par = para.cores, apply = lapply)
            if (verbose) 
                cat("NOT PARALLELIZING\n")
        }
        else {
            mc <- has.multi && !multicore:::isChild() && X > 
                1 && !is.na(para.cores) && (is.numeric(para.cores) && 
                para.cores > 1) || (is.logical(para.cores) && 
                para.cores == TRUE)
            par <- para.cores
            out.apply <- lapply
            if (mc) {
                if (is.logical(par) && par == TRUE) 
                  par <- multicore:::detectCores()
                par <- min(c(X, par, multicore:::detectCores()))
                if (verbose) 
                  cat("PARALLELIZING:", par, ": ")
                foreach.register.backend(par)
                if (verbose) 
                  cat(getDoParName(), getDoParWorkers(), "\n")
                out.apply <- function(list, FUN, ...) foreach(l = list) %dopar% 
                  {
                    FUN(l, ...)
                  }
            }
            else {
                par <- 1
                if (verbose) 
                  cat("NOT PARALLELIZING:", par, "\n")
            }
            out <- list(mc = mc, par = par, apply = out.apply)
        }
    }
    if (is.numeric(out$par) && !is.na(out$par)) 
        options(cores = out$par)
    else if (is.na(out$par) || (is.logical(out$par) && out$par == 
        TRUE)) 
        options(cores = NULL)
    else options(cores = 1)
    out
}
get.predictome.links <-
function (org.id = organism) 
{
    out <- list()
    for (i in c("chromo", "comp", "fusion", "phylogenetic")) {
        fname <- paste("data/predictome/predictome_", i, "_links.txt", 
            sep = "")
        pred.file <- paste("http://predictome.bu.edu/data/all", 
            i, "links.txt", sep = "_")
        err <- dlf(fname, pred.file, paste("Reading in predictome links from", 
            pred.file))
        pred.tab <- read.delim(gzfile(fname), head = T, as.is = T)
        pred.tab <- pred.tab[pred.tab$species == org.id, ]
        out[[i]] <- data.frame(protein1 = pred.tab$orf_id_1, 
            protein2 = pred.tab$orf_id_2, combined_score = 1)
    }
    closeAllConnections()
    out
}
get.prolinks.links <-
function (org.id = genome.info$org.id$V1[1]) 
{
    fname <- paste("data/", rsat.species, "/prolinks_", gsub(" ", 
        "_", org.id), ".txt", sep = "")
    org.file <- paste("http://mysql5.mbi.ucla.edu/public/Genomes/", 
        gsub(" ", "_", org.id), ".txt", sep = "")
    err <- dlf(fname, org.file, paste("Fetching PROLINKS links from", 
        org.file))
    prol.tab <- read.delim(gzfile(fname), head = T, as.is = T)
    fname <- paste("data/prolinks_GeneID_Genename.txt", sep = "")
    id.file <- "http://mysql5.mbi.ucla.edu/public/reference_files/GeneID_Genename.txt"
    err <- dlf(fname, id.file, paste("Fetching PROLINKS genename ref. file from", 
        id.file))
    id.tab <- read.delim(gzfile(fname), head = T, as.is = T)
    merged.tab <- merge(merge(prol.tab, id.tab, by.x = "gene_id_a", 
        by.y = "gene_id"), id.tab, by.x = "gene_id_b", by.y = "gene_id")
    out <- list()
    for (i in unique(merged.tab$method)) {
        out[[i]] <- merged.tab[merged.tab$method == i, c("name.x", 
            "name.y", "confidence")]
        colnames(out[[i]]) <- c("protein1", "protein2", "combined_score")
    }
    closeAllConnections()
    out
}
get.pssm.scans <-
function (motifs = "ALL", seqs = e$genome.info$all.upstream.seqs[[1]], 
    seq.type = names(e$meme.scores)[1], p.cutoff = "0.0001") 
{
    mast.cmd <- "./progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.0001 -ev 9999 -comp"
    mast.cmd <- sprintf("./progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt %s -ev 9999 -comp", 
        p.cutoff)
    ks <- 1:e$k.clust
    if (motifs != "ALL" && !is.null(motifs)) 
        ks <- as.integer(gsub("BIC_", "", unlist(get.biclusters(motif = motifs))))
    out <- do.call(rbind, mclapply(ks, function(k) {
        print(k)
        if (is.null(e$meme.scores[[seq.type]][[k]]) || e$meme.scores[[seq.type]][[k]] == 
            "" || is.null(e$meme.scores[[seq.type]][[k]]$meme.out)) 
            return(NULL)
        tmp1 <- e$cluster.meme.motif.lines(k, logodds = T)
        if (length(tmp1) < 10) 
            return(NULL)
        tmp3 <- e$runMast(tmp1, mast.cmd, names(seqs), seqs, 
            verbose = F, seq.type = seq.type, bg.list = NULL, 
            bg.fname = NULL, unlink = T)
        if (length(tmp3) <= 0) 
            return(NULL)
        tmp4 <- e$getMastPValuesAndEValues(tmp3, names(seqs))[[2]]
        if (nrow(tmp4) <= 0) 
            return(NULL)
        tmp4 <- cbind(tmp4, bic = rep(k, nrow(tmp4)))
        tmp4
    }))
    out
}
get.pv.ev.single <-
function (mast.out, rows) 
{
    pv.ev <- NULL
    if (length(grep("Error reading log-odds matrix file", mast.out)) <= 
        0 && class(mast.out) != "try-error" && length(mast.out) > 
        0) {
        pv.ev <- getMastPValuesAndEValues(mast.out, get.p.values = rows)
        attr(pv.ev, "mast.command.line") <- attr(mast.out, "mast.command.line")
        if (length(pv.ev) > 0 && nrow(pv.ev[[1]]) == 0 && nrow(pv.ev[[2]]) == 
            0) {
            pv.ev <- NULL
        }
        else {
            for (i in 1) {
                tmp <- as.matrix(pv.ev[[i]][, 2:ncol(pv.ev[[i]])])
                rownames(tmp) <- pv.ev[[i]][, 1]
                pv.ev[[i]] <- tmp
            }
        }
    }
    pv.ev
}
get.row.scores <-
function (k, cols = get.cols(k), for.rows = "all", ratios = ratios[[1]], 
    method = c("cor2", "abscor", "cor", "dist", "orig")[1], ...) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else rows <- k
    if (is.null(for.rows) || for.rows[1] == "all") 
        for.rows <- rownames(ratios)
    rows <- rows[rows %in% rownames(ratios)]
    cols <- cols[cols %in% colnames(ratios)]
    if (length(rows) < 1 || length(cols) < 1) 
        return(rep(NA, length(for.rows)))
    if (method == "orig") {
        rats <- ratios[for.rows, cols, drop = F]
        rats.mn <- colMeans(rats[rows, , drop = F], na.rm = T)
        rats[, ] <- t(t(rats) - rats.mn)^2
        col.weights <- if (exists("get.col.weights")) 
            get.col.weights(rows, cols, ratios)
        else NA
        if (is.na(col.weights[1])) 
            rats <- rowMeans(rats, na.rm = T)
        else rats <- apply(rats, 1, weighted.mean, w = col.weights[cols], 
            na.rm = T)
        rats <- log(rats + 1e-99)
    }
    else if (method == "pval") 
        rats <- get.row.scores.pVals(k, cols, rows, ratios, method, 
            ...)
    else if (exists("get.row.scores.NEW")) 
        rats <- get.row.scores.NEW(k, cols, rows, ratios, method, 
            ...)
    return(rats)
}
get.row.scores.NEW <-
function (k, cols = get.cols(k), rows, ratios = ratios[[1]], 
    method = c("cor2", "abscor", "cor", "dist")[1], ...) 
{
    rats <- ratios[, cols, drop = F]
    if (method == "dist" || substr(method, 1, 5) == "dist.") {
        require(proxy)
        if (method == "dist") 
            method <- "dist.Euclidean"
        d <- proxy::dist(rats, rats[rows, , drop = F], method = substring(method, 
            6))
        d[cbind(rows, rows)] <- NA
        rats <- log(apply(d, 1, mean, na.rm = T, ...) + 1e-99)
    }
    else if (method %in% c("cor", "cor2", "abscor")) {
        rats <- t(rats)
        d <- t(cor(rats[, rows, drop = F], rats, use = "pairwise", 
            ...))
        d[cbind(rows, rows)] <- NA
        if (method == "cor2") 
            d <- d^2
        else if (method == "abscor") 
            d <- abs(d)
        rats <- apply(log(1 - d + 1e-99), 1, mean, na.rm = T, 
            ...)
    }
    else if (method == "mi" || substr(method, 1, 3) == "mi.") {
        require(minet)
        if (method == "mi") 
            method <- "mi.spearman"
        d <- build.mim(t(rats), estimator = substring(method, 
            4), disc = "equalwidth")
        diag(d) <- NA
        rats <- -apply(d[, rows, drop = F], 1, mean, na.rm = T, 
            ...)
    }
    return(rats)
}
get.row.weights <-
function (rows, cols, ratios) 
NA
get.rows <-
function (k) 
clusterStack[[k]]$rows
get.sequence.psps <-
function (seqs) 
NULL
get.sequences <-
function (k, seq.type = paste(c("upstream", "upstream.noncod", 
    "upstream.noncod.same.strand", "downstream", "gene")[1], 
    "meme"), verbose = F, filter = T, distance = motif.upstream.search[[seq.type]], 
    ...) 
{
    if (length(k) <= 0) 
        return(NULL)
    if (is.numeric(k[1])) 
        rows <- get.rows(k)
    else if (!is.null(genome.info$genome.seqs) && k %in% names(genome.info$genome.seqs)) 
        return(genome.info$genome.seqs[k])
    else rows <- k
    if (is.null(rows)) 
        return(NULL)
    start.stops <- NULL
    if (is.na(seq.type) || strsplit(seq.type, " ")[[1]][1] == 
        "gene") 
        op.shift <- FALSE
    n.seq.type <- strsplit(seq.type, " ")[[1]][1]
    if (substr(n.seq.type, 1, 8) == "fstfile=") {
        if (!is.null(genome.info$all.upstream.seqs[[seq.type]]) && 
            length(genome.info$all.upstream.seqs[[seq.type]]) > 
                0) {
            seqs <- genome.info$all.upstream.seqs[[seq.type]]
        }
        else {
            seqs <- read.fasta(fname = gzfile(strsplit(n.seq.type, 
                "=")[[1]][2]))
        }
        seqs <- seqs[rows]
        names(seqs) <- toupper(rows)
    }
    else if (substr(n.seq.type, 1, 8) == "csvfile=") {
        if (!is.null(genome.info$all.upstream.seqs[[seq.type]]) && 
            length(genome.info$all.upstream.seqs[[seq.type]]) > 
                0) {
            seqs <- genome.info$all.upstream.seqs[[seq.type]]
        }
        else {
            tab <- read.csv(gzfile(strsplit(n.seq.type, "=")[[1]][2]), 
                head = F)
            seqs <- as.character(tab[, 2])
            names(seqs) <- toupper(as.character(tab[, 1]))
        }
        seqs <- seqs[rows]
        names(seqs) <- toupper(rows)
    }
    else {
        if (is.null(genome.info$feature.tab) || !"genome.seqs" %in% 
            names(genome.info) || is.null(genome.info$genome.seqs)) {
            if (!is.null(genome.info$all.upstream.seqs[[seq.type]])) {
                out.seqs <- genome.info$all.upstream.seqs[[seq.type]][rows]
                out.seqs <- out.seqs[!is.na(out.seqs) & out.seqs != 
                  ""]
                if (filter) 
                  out.seqs <- filter.sequences(out.seqs, NULL, 
                    seq.type, distance, verbose = verbose, ...)
                return(invisible(out.seqs))
            }
            else {
                stop("Motif searching is on but no ", seq.type, 
                  " sequences!")
            }
        }
        op.shift <- operon.shift[seq.type]
        coos <- get.gene.coords(rows, op.shift)
        if (is.null(coos) || nrow(coos) <= 0) 
            return(NULL)
        coos <- subset(coos, !is.na(start_pos) & !is.na(end_pos))
        if (is.null(coos) || nrow(coos) <= 0) 
            return(NULL)
        seqs <- character()
        if (n.seq.type %in% c("upstream.noncod", "upstream.noncod.same.strand")) {
            all.coos <- genome.info$feature.tab[, c("id", "name", 
                "contig", "strand", "start_pos", "end_pos")]
            all.coos <- subset(all.coos, name %in% unlist(genome.info$synonyms))
        }
        mc <- get.parallel(nrow(coos))
        tmp <- mc$apply(1:nrow(coos), function(i) {
            if (n.seq.type == "gene") {
                st.st <- coos[i, c("start_pos", "end_pos"), drop = F]
            }
            else if (n.seq.type == "upstream") {
                st.st <- if (coos$strand[i] == "D") 
                  c(coos$start_pos[i] - 1 - distance[2], coos$start_pos[i] - 
                    1 - distance[1])
                else c(coos$end_pos[i] + 1 + distance[1], coos$end_pos[i] + 
                  1 + distance[2])
            }
            else if (n.seq.type == "downstream") {
                st.st <- if (coos$strand[i] == "D") 
                  c(coos$end_pos[i] + 1 + distance[1], coos$end_pos[i] + 
                    1 + distance[2])
                else c(coos$start_pos[i] - 1 - distance[2], coos$start_pos[i] - 
                  1 - distance[1])
            }
            else if (n.seq.type %in% c("upstream.noncod", "upstream.noncod.same.strand")) {
                cc <- all.coos[as.character(all.coos$contig) == 
                  as.character(coos$contig[i]) & abs(all.coos$start_pos - 
                  coos$start_pos[i]) <= 1e+05, ]
                if (n.seq.type == "upstream.noncod.same.strand") 
                  cc <- all.coos[as.character(all.coos$strand) == 
                    as.character(coos$strand[i]), ]
                if (coos$strand[i] == "D") {
                  nearest <- max(cc$end_pos[cc$end_pos < coos$start_pos[i]])
                  st.st <- c(nearest, coos$start_pos[i] - distance[1] - 
                    1)
                }
                else if (coos$strand[i] == "R") {
                  nearest <- min(cc$start_pos[cc$start_pos > 
                    coos$end_pos[i]])
                  st.st <- c(coos$end_pos[i] + distance[1] + 
                    1, nearest)
                }
            }
            seq <- substr(genome.info$genome.seqs[[as.character(coos$contig[i])]], 
                st.st[1], st.st[2])
            if (coos$strand[i] == "R") 
                seq <- rev.comp(seq)
            if (nchar(seq) > abs(diff(distance))) {
                if (coos$strand[i] == "D") 
                  seq <- substr(seq, 1, abs(diff(distance)))
                else seq <- rev.comp(substr(rev.comp(seq), 1, 
                  abs(diff(distance))))
            }
            out <- list(seq = seq, name = as.character(coos$names[i]), 
                start.stops = data.frame(start = st.st[1], end = st.st[2], 
                  strand = as.character(coos$strand[i]), contig = as.character(coos$contig[i])))
            out
        })
        for (i in tmp) {
            seqs[i$name] <- i$seq
            start.stops <- rbind(start.stops, i$start.stops)
            rownames(start.stops)[nrow(start.stops)] <- i$name
        }
        rownames(start.stops) <- names(seqs) <- make.unique(rownames(start.stops))
        rows <- rows[rows %in% names(seqs)]
        start.stops <- start.stops[rows, , drop = F]
        seqs <- seqs[rows]
        names(seqs) <- rownames(start.stops) <- rows
    }
    if (any(is.na(seqs))) {
        warning("Warning: could not find '", n.seq.type, "' sequences for all input genes", 
            call. = F)
        if (!is.null(start.stops)) 
            start.stops <- start.stops[!is.na(seqs), ]
        seqs <- seqs[!is.na(seqs)]
    }
    if (filter) 
        seqs <- filter.sequences(seqs, start.stops, seq.type, 
            distance, verbose = verbose, ...)
    attr(seqs, "start.stops") <- start.stops
    invisible(seqs)
}
get.stats <-
function (mean.func = mean) 
{
    changed <- NA
    if (!exists("row.memb")) {
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        if (is.vector(row.memb)) 
            row.memb <- t(row.memb)
        rownames(row.memb) <- attr(ratios, "rnames")
        col.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "cnames") %in% get.cols(k))
        if (is.vector(col.memb)) 
            col.memb <- t(col.memb)
        rownames(col.memb) <- attr(ratios, "cnames")
    }
    if (!exists("row.scores") || is.null(row.scores)) {
        if (attr(get.all.scores, "version") == 1) {
            tmp <- get.all.scores()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
        }
        else if (attr(get.all.scores, "version") == 2) {
            tmp <- get.old.scores.matrices()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
        }
    }
    cs <- as.list(clusterStack)
    resids <- sapply(cs, "[[", "resid")
    if (is.matrix(resids)) 
        resids <- apply(resids, 1, function(r) mean.func(r[r != 
            1], na.rm = T))
    else resids <- mean.func(resids[resids != 1], na.rm = T)
    p.clusts <- sapply(cs, "[[", "p.clust")
    if (is.matrix(p.clusts)) 
        p.clusts <- apply(p.clusts, 1, mean.func, na.rm = T)
    else p.clusts <- mean.func(p.clusts, na.rm = T)
    out <- data.frame(iter = iter, changed = changed, row.scores = mean.func(row.scores[, 
        ][row.memb[, ]], na.rm = T), col.scores = mean.func(col.scores[, 
        ][col.memb[, ]], na.rm = T), mot.scores = if (!is.null(mot.scores)) 
        mean.func(mot.scores[, ][row.memb[, ]], na.rm = T)
    else NA, net.scores = if (!is.null(net.scores)) 
        mean.func(net.scores[, ][row.memb[, ]], na.rm = T)
    else NA, resid = weighted.mean(resids, row.weights, na.rm = T), 
        nrow = mean.func(sapply(cs, "[[", "nrows"), na.rm = T), 
        ncol = mean.func(sapply(cs, "[[", "ncols"), na.rm = T), 
        p.clust = if (!all(is.na(p.clusts))) 
            weighted.mean(p.clusts, mot.weights, na.rm = T)
        else NA)
    if (length(resids) > 1) 
        for (i in names(resids)) {
            out <- cbind(out, resids[i])
            names(out)[ncol(out)] <- paste("resid", i, sep = ".")
        }
    if (length(p.clusts) > 1) 
        for (i in names(p.clusts)) {
            out <- cbind(out, p.clusts[i])
            names(out)[ncol(out)] <- paste("p.clust", i, sep = ".")
        }
    if (length(networks) > 1) {
        for (i in names(net.weights)) {
            if (exists("cluster.net.scores") && i %in% colnames(cluster.net.scores)) 
                out <- cbind(out, weighted.mean(cluster.net.scores[, 
                  i], sapply(cs, "[[", "nrows"), na.rm = T))
            else out <- cbind(out, rep(NA, nrow(out)))
            names(out)[ncol(out)] <- paste("net", i, sep = ".")
        }
        if (exists("cluster.net.scores") && "net.scores" %in% 
            colnames(cluster.net.scores)) {
            out[, "net.scores"] <- weighted.mean(cluster.net.scores[, 
                "net.scores"], sapply(cs, "[[", "nrows"), na.rm = T)
        }
        else {
            out[, "net.scores"] <- mean.func(net.scores[, ][row.memb[, 
                ]], na.rm = T)
        }
    }
    out
}
get.synonyms <-
function (gns, ft = genome.info$feature.names, ignore.case = T, 
    verbose = F, fast = F, force = F) 
{
    if (exists("no.genome.info") && no.genome.info) {
        out <- as.list(gns)
        names(out) <- gns
        return(out)
    }
    out <- list()
    if ((!force && exists("genome.info") && !is.null(genome.info$synonyms))) {
        gns.cached <- gns[gns %in% names(genome.info$synonyms)]
        out <- genome.info$synonyms[gns.cached]
        gns <- gns[!gns %in% names(genome.info$synonyms)]
        if (length(gns) <= 0 || (is.null(ft) && (!exists("translation.tab") || 
            is.null(translation.tab)))) 
            return(out)
    }
    tmp.out <- as.list(gns)
    names(tmp.out) <- gns
    if (is.null(ft) && (!exists("translation.tab") || is.null(translation.tab))) 
        return(c(out, tmp.out))
    gns.orig <- gns
    gns <- gsub("m$|_\\d$|\\-\\S$", "", gns, perl = T)
    gns <- gsub("([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", "\\\\\\1", 
        gns, perl = T)
    gns <- gns[!is.na(gns) & gns != ""]
    ft <- ft[, c("id", "names")]
    if (exists("translation.tab") && !is.null(translation.tab)) 
        ft <- rbind(ft, data.frame(id = as.character(translation.tab$V1), 
            names = as.character(translation.tab$V2)))
    ft <- subset(ft, names != "")
    if (verbose) 
        ggggg <- gns[seq(100, length(gns), by = 100)]
    mc <- get.parallel(length(gns), verbose = F)
    tmp <- mc$apply(gns, function(g) {
        if (verbose && g %in% ggggg) 
            cat(" ...", g)
        greg <- paste("^", g, sep = "")
        tmp <- subset(ft, grepl(greg, id, perl = T, ignore = ignore.case) | 
            grepl(greg, names, perl = T, ignore = ignore.case))
        if (nrow(tmp) <= 0) 
            return(g)
        tmp2 <- unique(c(g, as.character(tmp$id), as.character(tmp$names)))
        if (!fast) {
            tmp2 <- subset(ft, id %in% tmp2 | names %in% tmp2)
            tmp2 <- unique(c(g, as.character(tmp2[, 1]), as.character(tmp2[, 
                2])))
        }
        tmp2 <- gsub("\\\\([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", 
            "\\1", tmp2, perl = T)
        tmp2
    })
    names(tmp) <- gns.orig
    if (verbose) 
        cat("\n")
    c(tmp, out)
}
get.type <-
function (obj) 
{
    prefix <- sapply(strsplit(obj, "_"), "[", 1)
    out <- ifelse(prefix == "BIC", "bicluster", ifelse(prefix == 
        "MOT", "motif", ifelse(prefix == "MOTC", "motif.cluster", 
        ifelse(obj %in% attr(e$ratios, "rnames"), "gene", "condition"))))
    out
}
get.unpreprocessed.ratios <-
function (...) 
{
    return(ratios.raw)
}
get.updated.memberships <-
function (row.membership, col.membership, rr.scores, cc.scores) 
{
    rm <- t(apply(rr.scores, 1, order, decreasing = T)[1:n.clust.per.row, 
        , drop = F])
    rm <- t(apply(rm, 1, sort))
    if (n.clust.per.row == 1) 
        rm <- t(rm)
    if (ncol(rm) < ncol(row.membership)) 
        rm <- cbind(rm, matrix(0, nrow = nrow(rm), ncol = ncol(row.membership) - 
            ncol(rm)))
    for (i in 1:nrow(rm)) {
        if (all(rm[i, ] %in% row.membership[i, ])) 
            next
        mc <- max.changes["rows"]
        if (mc < 1 && mc > 0 && runif(1) > mc) 
            next
        for (ii in 1:mc) {
            if (sum(!rm[i, ] %in% row.membership[i, ]) >= mc) {
                if (any(row.membership[i, ] == 0)) {
                  col.change <- which(row.membership[i, ] == 
                    0)[1]
                }
                else {
                  ttmp <- tabulate(row.membership[i, ])
                  if (any(ttmp > 1)) {
                    col.change <- which(row.membership[i, ] %in% 
                      which(ttmp > 1))[1]
                  }
                  else {
                    delta <- rr.scores[i, rm[i, ]] - rr.scores[i, 
                      row.membership[i, ]]
                    if (any(row.membership[i, ] %in% rm[i, ])) 
                      delta[row.membership[i, ] %in% rm[i, ]] <- 0
                    if (all(is.na(delta) | delta <= 0)) 
                      next
                    col.change <- which.max(delta)
                  }
                }
                if (exists("maintain.seed") && !is.null(maintain.seed) && 
                  !is.null(maintain.seed$rows) && !is.null(maintain.seed$rows[[as.character(row.membership[i, 
                  col.change])]]) && rownames(row.membership)[i] %in% 
                  maintain.seed$rows[[as.character(row.membership[i, 
                    col.change])]]) 
                  next
                if (!rm[i, col.change] %in% row.membership[i, 
                  ]) 
                  row.membership[i, col.change] <- rm[i, col.change]
            }
        }
    }
    if (!is.null(cc.scores)) {
        cm <- t(apply(cc.scores, 1, order, decreasing = T)[1:n.clust.per.col, 
            , drop = F])
        if (ncol(cm) < ncol(col.membership)) 
            cm <- cbind(cm, matrix(0, nrow = nrow(cm), ncol = ncol(col.membership) - 
                ncol(cm)))
        for (i in 1:nrow(cm)) {
            mc <- max.changes["cols"]
            if (mc < 1 && mc > 0 && runif(1) > mc) 
                next
            for (ii in 1:mc) {
                if (sum(!cm[i, ] %in% col.membership[i, ]) >= 
                  mc) {
                  if (any(col.membership[i, ] == 0)) {
                    col.change <- which(col.membership[i, ] == 
                      0)[1]
                  }
                  else {
                    ttmp <- tabulate(col.membership[i, ])
                    if (any(ttmp > 1)) {
                      col.change <- which(col.membership[i, ] %in% 
                        which(ttmp > 1))[1]
                    }
                    else {
                      delta <- cc.scores[i, cm[i, ]] - cc.scores[i, 
                        col.membership[i, ]]
                      if (all(is.na(delta) | delta <= 0)) 
                        next
                      col.change <- which.max(delta)
                    }
                  }
                  if (exists("maintain.seed") && !is.null(maintain.seed) && 
                    !is.null(maintain.seed$cols) && !is.null(maintain.seed$cols[[as.character(col.membership[i, 
                    col.change])]]) && colnames(col.membership)[i] %in% 
                    maintain.seed$cols[[as.character(col.membership[i, 
                      col.change])]]) 
                    next
                  col.membership[i, col.change] <- cm[i, col.change]
                }
            }
        }
    }
    invisible(list(r = row.membership, c = col.membership))
}
get.updated.memberships.NEWTRY <-
function (rows = TRUE) 
{
    tmp2 = which(r.scores < Inf, arr = T)
    tmp2 = cbind(tmp2, r.scores[tmp2])
    tmp2 = tmp2[order(tmp2[, 3]), ]
    nr <- nrow(row.membership)
    rm <- matrix(0, nrow = nr, ncol = max(table(tmp2[1:(nr * 
        n.clust.per.row), 1])))
    rownames(rm) <- rownames(row.membership)
    for (i in 1:(attr(ratios, "nrow") * n.clust.per.row)) {
        r <- rownames(tmp2)[i]
        rm[r, min(which(rm[r, ] == 0))] <- tmp2[i, 2]
    }
    tmp2 = which(c.scores < Inf, arr = T)
    tmp2 = cbind(tmp2, c.scores[tmp2])
    tmp2 = tmp2[order(tmp2[, 3]), ]
    nr <- nrow(col.membership)
    cm <- matrix(0, nrow = nr, ncol = max(table(tmp2[1:(nr * 
        n.clust.per.col), 1])))
    rownames(cm) <- rownames(col.membership)
    for (i in 1:(attr(ratios, "ncol") * n.clust.per.col)) {
        r <- rownames(tmp2)[i]
        cm[r, min(which(cm[r, ] == 0))] <- tmp2[i, 2]
    }
    cm <- cbind(cm, col.membership)
    cm <- apply(cm, 1, function(i) c(unique(i), rep(0, ncol(cm) - 
        length(unique(i)))))
    invisible(list(r = rm, c = cm))
}
getMastPValuesAndEValues <-
function (mastOutput, get.p.values = NULL) 
{
    lines <- grep("COMBINED P-VALUE", mastOutput)
    if (length(lines) > 0) {
        splitted <- strsplit(mastOutput[lines], "[\\t\\s]+", 
            perl = T)
        out <- t(sapply(1:length(lines), function(i) {
            gene <- mastOutput[lines[i] - 2]
            splt <- splitted[[i]]
            p.val <- splt[8]
            e.val <- splt[11]
            c(gene = gene, p.value = p.val, e.value = e.val)
        }))
        out <- data.frame(gene = out[, "gene"], p.value = as.numeric(out[, 
            "p.value"]), e.value = as.numeric(out[, "e.value"]))
    }
    out2 <- data.frame()
    if (!is.null(get.p.values) && !is.na(get.p.values)) {
        tmp <- get.mast.pvals(mastOutput, in.genes = get.p.values)
        for (g in names(tmp)) {
            pv <- as.numeric(tmp[[g]]$pvals)
            pos <- as.integer(tmp[[g]]$posns)
            mots <- as.integer(tmp[[g]]$mots)
            if (!all(c(length(pv), length(pos)) == length(mots))) 
                pv <- c(pv, rep(pv[1], length(pos) - length(pv)))
            out2 <- rbind(out2, data.frame(gene = g, pvals = pv, 
                posns = pos, mots = mots))
        }
    }
    return(list(out, out2))
}
getMemeMotifInfo <-
function (memeOutput) 
{
    out <- list()
    lines <- grep("^MOTIF\\s+\\d", memeOutput, perl = T)
    if (length(lines) <= 0) 
        lines <- grep("^MOTIF\\s+", memeOutput, perl = T)
    if (length(lines) > 0) {
        pssms <- getMemeMotifPssm(memeOutput, n.motif = length(lines))
        splitted <- strsplit(memeOutput[lines], "[\\t\\s]+", 
            perl = T)
        for (i in 1:length(lines)) {
            splt <- splitted[[i]]
            motif <- as.integer(splt[2])
            width <- as.integer(splt[5])
            sites <- as.integer(splt[8])
            llr <- as.integer(splt[11])
            e.value <- as.numeric(sub("\\+", "", splt[14]))
            pssm <- pssms[[motif]]$pssm
            l2 <- grep(paste("Motif", motif, "sites sorted by position p-value"), 
                memeOutput) + 4
            l3 <- grep("--------------------------------------------------------------------------------", 
                memeOutput[(l2 + 1):length(memeOutput)])[1] + 
                l2 - 1
            posns <- do.call(rbind, strsplit(memeOutput[l2:l3], 
                "[\\t\\s]+", perl = T))[, c(1:4, 6)]
            colnames(posns) <- c("gene", "strand", "start", "p.value", 
                "site")
            posns <- data.frame(gene = posns[, "gene"], strand = posns[, 
                "strand"], start = as.integer(posns[, "start"]), 
                p.value = as.numeric(posns[, "p.value"]), site = posns[, 
                  "site"])
            out[[motif]] <- list(width = width, sites = sites, 
                llr = llr, e.value = e.value, pssm = pssm, posns = posns)
        }
    }
    out
}
getMemeMotifPssm <-
function (memeOut, n.motif = 1) 
{
    pssms <- list()
    for (i in 1:n.motif) {
        m.line1 <- grep(sprintf("Motif %d position-specific probability matrix", 
            i), memeOut)
        if (length(m.line1) > 0) {
            m.desc <- strsplit(memeOut[m.line1 + 2], " ")[[1]]
            winLen <- as.numeric(m.desc[6])
            e.val <- as.numeric(m.desc[10])
            pssm <- do.call(rbind, strsplit(memeOut[m.line1 + 
                2 + 1:winLen], "\\s+", perl = T))[, 2:5]
            pssm <- matrix(as.numeric(pssm), nrow = winLen, ncol = 4, 
                byrow = F)
            pssms[[i]] <- list(pssm = pssm, e.val = e.val)
        }
        else {
            pssms[[i]] <- list(pssm = NULL, e.val = 99999)
        }
    }
    return(pssms)
}
gibbs.one.cluster <-
function (k, seq.type = names(mot.weights)[1], n.motifs = 2, 
    width = 8, verbose = F, unlink = T, ...) 
{
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    seqs <- get.sequences(rows, seq.type = seq.type, ...)
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    bgs <- genome.info$all.upstream.seqs[[seq.type]]
    fname <- my.tempfile("gibbs_")
    cat(paste(">", names(seqs), "\n", seqs, sep = ""), file = fname, 
        sep = "\n")
    bgtab <- data.frame(names(genome.info$bg.list[[seq.type]]), 
        sapply(genome.info$bg.list[[seq.type]], "[", 1))
    rownames(bgtab) <- NULL
    colnames(bgtab) <- c("V1", "V2")
    bgtab <- bgtab[-1, ]
    bgtab$V1 <- tolower(bgtab$V1)
    bgfname <- my.tempfile("gibbs_bg_")
    write.table(bgtab, bgfname, sep = "\t", quote = F, col.names = F, 
        row.names = F)
    outfname <- my.tempfile("gibbs_out_")
    cmd <- sprintf("%s/gibbs/Gibbs.i686-apple-darwin %s %d %d -n", 
        progs.dir, fname, width, floor(sqrt(length(seqs))))
    if (verbose) 
        print(cmd)
    out <- system(cmd, intern = T, ignore.stderr = !verbose)
    out <- readLines(outfname)
    if (unlink) 
        unlink(c(fname, bgfname, outfname))
    out
}
glam.one.cluster <-
function (k, seq.type = "upstream", min.seqs = cluster.rows.allowed[1], 
    max.seqs = cluster.rows.allowed[2], verbose = T, unlink = T, 
    ...) 
{
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    seqs <- get.sequences(rows, seq.type = seq.type, ...)
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k))
    all.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    tmp <- strsplit(meme.cmd, " ")[[1]]
    w.min <- tmp[which(tmp == "-minw") + 1]
    w.max <- tmp[which(tmp == "-maxw") + 1]
    fna.file <- my.tempfile("glamseqs.fna.")
    fna1.file <- my.tempfile("glamseqs1.fna.")
    all.fna.file <- my.tempfile("glamseqs.all.fna.")
    out.file <- my.tempfile("glam.out.")
    glam.out <- list()
    cat(paste(">", names(seqs), "\n", seqs, sep = ""), file = fna.file, 
        sep = "\n", append = F)
    file.copy(fna.file, fna1.file, overwrite = T)
    cat(paste(">", names(all.seqs), "\n", all.seqs, sep = ""), 
        file = all.fna.file, sep = "\n", append = F)
    for (i in 1:n.motifs[[seq.type]][iter]) {
        cmd <- paste(sprintf("%s/glam2 -2 -n 5000 -a", progs.dir), 
            w.min, "-b", w.max, "-q 1 -O", out.file, "n", fna1.file)
        if (verbose) 
            cat(i, cmd, "\n")
        system(cmd)
        glam.out[[i]] <- list()
        glam.out[[i]]$txt <- readLines(paste(out.file, "/glam2.txt", 
            sep = ""))
        glam.out[[i]]$meme <- readLines(paste(out.file, "/glam2.meme", 
            sep = ""))
        cmd <- paste(sprintf("%s/glam2scan -n 5000 -2 n ", progs.dir), 
            out.file, "/glam2.txt ", all.fna.file, sep = "")
        if (verbose) 
            cat(i, cmd, "\n")
        glam.out[[i]]$scan <- system(cmd, intern = T)
        if (i < n.motifs[[seq.type]][iter]) {
            cmd <- paste(sprintf("%s/glam2mask -o ", progs.dir), 
                fna1.file, " ", out.file, "/glam2.txt ", fna.file, 
                sep = "")
            if (verbose) 
                cat(i, cmd, "\n")
            system(cmd)
            file.copy(fna1.file, fna.file, overwrite = T)
        }
    }
    if (unlink) 
        system(paste("rm -rf", fna.file, fna1.file, all.fna.file, 
            out.file))
    for (i in 1:length(glam.out)) {
        tmp <- glam.out[[i]]$scan
        tmp <- tmp[tmp != ""]
        tmp <- tmp[!1:length(tmp) %in% grep("^\\s+", tmp, perl = T)]
        tmp <- tmp[-(1:3)]
        tmp <- do.call(rbind, strsplit(tmp, "\\s+"))
        scores <- as.numeric(tmp[, 6])
        names(scores) <- tmp[, 1]
        glam.out[[i]]$posns <- data.frame(gene = tmp[, 1], start = as.integer(tmp[, 
            2]), end = as.integer(tmp[, 4]), strand = tmp[, 5], 
            score = scores)
    }
    for (i in 1:length(glam.out)) {
        tmp <- glam.out[[i]]$txt
        score <- as.numeric(strsplit(grep("^Score:\\s+", tmp, 
            perl = T, val = 1)[1], "\\s+", perl = T)[[1]][2])
        start <- grep("a  c  g  t Del Ins Score", tmp, fixed = T)[1] + 
            1
        end <- grep("^Score:\\s+", tmp, perl = T)[2] - 2
        tmp <- tmp[start:end]
        tmp <- gsub("^\\s+", "", tmp[seq(1, length(tmp), by = 2)], 
            perl = T)
        tmp <- do.call(rbind, strsplit(tmp, "\\s+"))[, 1:4]
        tmp <- t(apply(tmp, 1, as.numeric))
        end <- grep("^Score:\\s+", tmp, perl = T)[2] - 2
        attr(tmp, "score") <- score
        colnames(tmp) <- c("A", "C", "G", "T")
        glam.out[[i]]$pssm <- tmp
    }
    for (i in 1:length(glam.out)) {
        attr(glam.out[[i]], "meme.out") <- glam.out[[i]]$meme
        glam.out[[i]]$txt <- glam.out[[i]]$scan <- glam.out[[i]]$meme <- NULL
    }
    invisible(glam.out)
}
id.duplicate.clusters <-
function (scores = r.scores, cor.cutoff = 0.9) 
{
    cors <- cor(scores[, ], use = "pairwise", method = "pearson")
    cors[lower.tri(cors, diag = T)] <- NA
    tmp <- which(cors >= cor.cutoff, arr = T)
    cbind(tmp, cors[tmp])
}
install.binaries <-
function (meme.version = "4.3.0", url = sprintf("http://meme.nbcr.net/downloads/old_versions/meme_%s.tar.gz", 
    meme.version), make = "make -j 4") 
{
    cwd <- setwd(system.file(package = "cMonkey"))
    on.exit(setwd(cwd))
    if (!exists("progs")) 
        dir.create("progs")
    setwd("progs/")
    cMonkey:::dlf(sprintf("meme_%s.tar.gz", meme.version), url)
    system(sprintf("tar -xzf meme_%s.tar.gz", meme.version))
    unlink(sprintf("meme_%s.tar.gz", meme.version))
    setwd(sprintf("meme_%s", meme.version))
    dir.create("local")
    system(sprintf("./configure --prefix=%s/local/ --enable-dependency-tracking --enable-opt --disable-shared --disable-fast-install --enable-serial --enable-build-libxml2 --enable-build-libxslt --disable-shared --enable-static --with-gnu-ld", 
        getwd()))
    system(make)
    system("make install")
    setwd("..")
    system(sprintf("ln -s meme_%s/local/bin/meme", meme.version))
    system(sprintf("ln -s meme_%s/local/bin/mast", meme.version))
    system(sprintf("ln -s meme_%s/local/bin/dust", meme.version))
    system(sprintf("ln -s meme_%s/local/bin/tomtom", meme.version))
    setwd(cwd)
}
list.reference <-
function (l, file, ...) 
{
    if (!big.memory || !require(filehash)) 
        return(l)
    if (!is.null(l) && ("filehashDB1" %in% class(l) || length(l) <= 
        0)) 
        return(l)
    if (big.memory && !file.exists(cmonkey.filename)) 
        dir.create(cmonkey.filename, recursive = T, show = F)
    if (big.memory.verbose) 
        try(message("Filehashing: ", file), silent = T)
    if (file.exists(file)) 
        unlink(file, recursive = T)
    if (is.null(l)) {
        dbCreate(file, ...)
        l <- dbInit(file, ...)
    }
    else l <- dumpList(l, file, ...)
    l
}
load.genome.info.MicrobesOnline <-
function () 
{
    f <- sprintf("data/%s/microbesOnlineGenomeInfo_%d.tsv", rsat.species, 
        taxon.id)
    if (!file.exists(sprintf("%s.gz", f))) {
        try(dlf(f, sprintf("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab", 
            taxon.id)))
        system(sprintf("gzip -fv %s", f))
    }
    cat("Loading:", f, "\n")
    out <- read.delim(gzfile(sprintf("%s.gz", f)), row.names = 1, 
        sep = "\t", as.is = T, head = T)
    out
}
load.ratios <-
function (ratios) 
{
    if (is.null(ratios)) 
        return(NULL)
    if (is.character(ratios) && file.exists(ratios)) {
        cat("Loading ratios file", ratios, "\n")
        ratios <- read.delim(file = gzfile(ratios), sep = "\t", 
            as.is = T, header = T)
    }
    if (is.matrix(ratios) || is.data.frame(ratios)) {
        if (class(ratios[, 1]) == "character") {
            ratios <- ratios[!duplicated(ratios[, 1]), ]
            rownames(ratios) <- attr(ratios, "rnames") <- ratios[, 
                1]
            ratios <- ratios[, -1]
        }
        if (class(ratios[, 1]) == "character") 
            ratios <- ratios[, -1]
    }
    cat("Original ratios matrix is", paste(dim(ratios), collapse = "x"), 
        "\n")
    if (!is.matrix(ratios)) 
        ratios <- as.matrix(ratios)
    if (is.null(attr(ratios, "isPreProcessed")) || attr(ratios, 
        "isPreProcessed") == FALSE) {
        ratios <- preprocess.ratios(ratios)
        attr(ratios, "isPreProcessed") <- TRUE
    }
    closeAllConnections()
    ratios
}
load.ratios.GEO <-
function (search.terms) 
{
    if (missing(search.terms)) 
        search.terms <- gsub("_", " ", rsat.species)
    message("Searching GEO for terms='", search.terms, "'")
    search.terms <- URLencode(search.terms)
    try(dlf(sprintf("data/%s/GEO_search.xml", rsat.species), 
        sprintf("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=%s&retmax=99999", 
            search.terms)))
    tmp <- readLines(sprintf("data/%s/GEO_search.xml", rsat.species))
    tmp <- grep("<Id>\\d+</Id>", tmp, perl = T, val = T)
    ids <- sapply(strsplit(tmp, "[<>]"), "[", 3)
    ids <- gsub("^20+", "", gsub("^10+", "", ids))
    cat("Downloading GEO series:", ids, "\n")
    for (id in ids) {
        err <- try(dlf(sprintf("data/%s/GSE%s_series_matrix.txt.gz", 
            rsat.species, id), sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE%s/GSE%s_series_matrix.txt.gz", 
            id, id)))
        if (class(err) == "try-error") {
            require(RCurl)
            tmp <- try(getURL(sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE%s/", 
                id)))
            i <- 0
            while (class(tmp) == "try-error" && (i <- i + 1) < 
                10) tmp <- try(getURL(sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE%s/", 
                id)))
            tmp <- strsplit(tmp, "\n")[[1]]
            tmp <- sapply(tmp, function(i) strsplit(i, "\\s+")[[1]][9])
            names(tmp) <- NULL
            for (id2 in tmp) {
                err <- try(dlf(sprintf("data/%s/%s", rsat.species, 
                  id2), sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE%s/%s", 
                  id, id2)))
            }
        }
    }
    out <- list()
    for (f in list.files(patt = glob2rx("GSE*_series_matrix.txt.gz"), 
        path = sprintf("data/%s/", rsat.species), full = T)) {
        if (file.info(f)$size == 0) 
            next
        cat("Loading:", f, "\n")
        rats <- read.delim(gzfile(f), comment = "!")
        rownames(rats) <- sapply(strsplit(as.character(rats[[1]]), 
            "|", fixed = T), "[", 1)
        out[[basename(f)]] <- as.matrix(rats[, -1])
    }
    out
}
load.ratios.MicrobesOnline <-
function () 
{
    f <- sprintf("data/%s/microbesOnlineExprData_%d.tsv", rsat.species, 
        taxon.id)
    if (!file.exists(sprintf("%s.gz", f))) {
        timeout <- options(timeout = 600)
        try(dlf(f, sprintf("http://www.microbesonline.org/cgi-bin/microarray/getData.cgi?taxId=%d", 
            taxon.id)))
        system(sprintf("gzip -fv %s", f))
        options(timeout = timeout$timeout)
    }
    cat("Loading:", f, "\n")
    rats <- read.delim(gzfile(sprintf("%s.gz", f)), row.names = 1, 
        sep = "\t", as.is = T, head = T)
    as.matrix(rats)
}
load.sif.interactions <-
function (sif.fname) 
{
    sif <- read.delim(gzfile(sif.fname), sep = "", head = F, 
        comment = "#")
    contains.weights <- ncol(sif) == 3 && any(sapply(1:ncol(sif), 
        function(i) is.numeric(sif[, i])))
    if (contains.weights) {
        weight.col <- which(sapply(1:ncol(sif), function(i) is.numeric(sif[, 
            i])))
        if (length(weight.col) == 1) {
            sif[, weight.col] <- sif[, weight.col] * 1000/max(sif[, 
                weight.col], na.rm = T)
            if (weight.col == 1) 
                sif <- sif[, c(2, 1, 3)]
            else if (weight.col == 3) 
                sif <- sif[, c(1, 3, 2)]
            colnames(sif) <- c("V1", "V2", "V3")
        }
    }
    else if (ncol(sif) == 3) {
        sif <- data.frame(V1 = sif$V1, V2 = rep(1000, nrow(sif)), 
            V3 = sif$V3)
    }
    else {
        sif <- data.frame(V1 = sif$V1, V2 = rep(1000, nrow(sif)), 
            V3 = sif$V2)
    }
    sif$V2[is.na(sif$V2)] <- 0
    sif <- sif[, c("V1", "V3", "V2")]
    colnames(sif) <- c("protein1", "protein2", "combined_score")
    closeAllConnections()
    sif
}
make.pv.ev.matrix <-
function (out.ms, make.ev = F) 
{
    mot.rows <- character()
    for (k in 1:k.clust) {
        if (is.null(out.ms[[k]]$pv.ev)) 
            next
        mot.rows <- unique(c(mot.rows, rownames(out.ms[[k]]$pv.ev[[1]])))
    }
    mot.rows <- sort(mot.rows)
    out.pv <- out.ev <- NULL
    for (k in 1:k.clust) {
        m <- out.ms[[k]]
        if (is.null(m) || is.null(m$pv.ev)) {
            out.pv <- cbind(out.pv, rep(NA, length(mot.rows)))
            if (make.ev) 
                out.ev <- cbind(out.ev, rep(NA, length(mot.rows)))
        }
        else {
            m.scores <- numeric(length = length(mot.rows))
            tmp <- m$pv.ev[[1]][, "p.value"]
            names(tmp) <- rownames(m$pv.ev[[1]])
            m.scores <- tmp[mot.rows]
            out.pv <- cbind(out.pv, m.scores)
            colnames(out.pv) <- NULL
            if (make.ev) {
                m.scores <- numeric(length = length(mot.rows))
                tmp <- m$pv.ev[[1]][, "e.value"]
                names(tmp) <- rownames(m$pv.ev[[1]])
                m.scores <- tmp[mot.rows]
                out.ev <- cbind(out.ev, m.scores)
                colnames(out.ev) <- NULL
            }
            out.ms[[k]]$pv.ev[[1]] <- NULL
        }
    }
    rownames(out.pv) <- mot.rows
    if (!is.null(out.pv)) 
        rownames(out.pv) <- mot.rows
    out.pv
}
matrix.reference <-
function (m, ...) 
{
    if (!big.memory) {
        return(m)
    }
    if (!require(bigmemory) || is.big.matrix(m)) 
        return(m)
    options(bigmemory.allow.dimnames = TRUE)
    if (!file.exists(cmonkey.filename)) 
        dir.create(cmonkey.filename, recursive = T, show = F)
    if (!"backingfile" %in% names(list(...))) 
        backingfile <- as.character(substitute(m))
    else backingfile <- list(...)$backingfile
    if (!"backingpath" %in% names(list(...))) 
        backingpath <- cmonkey.filename
    else backingpath <- list(...)$backingpath
    if (big.memory.verbose) 
        message("Filebacking: ", backingpath, "/", backingfile)
    m <- as.big.matrix(m, backingfile = backingfile, backingpath = backingpath)
    m
}
meme.one.cluster <-
function (k, seq.type = names(mot.weights)[1], verbose = F, force = F, 
    ...) 
{
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    meme.out <- mast.out <- NULL
    cmd <- sprintf(meme.cmd[seq.type], n.motifs[[seq.type]][iter])
    pal.opt <- "non"
    if ("pal.opt" %in% names(list(...))) 
        pal.opt <- list(...)$pal.opt
    get.meme.consensus <- function(cmd, min.iter = 500, max.e.value = 0.1, 
        ...) {
        if (grepl("-cons $compute", cmd, fixed = T)) {
            if (iter > min.iter && !is.null(ms) && !is.null(ms$meme.out)) {
                e.val <- sapply(1:length(ms$meme.out), function(i) ms$meme.out[[i]]$e.value)
                if (min(e.val, na.rm = T) < max.e.value) {
                  best <- which.min(e.val)
                  consensus <- toupper(pssm.to.string(ms$meme.out[[best]]$pssm))
                  cmd <- gsub("$compute", consensus, cmd, fixed = T)
                }
            }
        }
        if (grepl("-cons $compute", cmd, fixed = T)) 
            cmd <- gsub("-cons $compute", "", cmd, fixed = T)
        cmd
    }
    ms <- try(meme.scores[[seq.type]][[k]])
    if ("try-error" %in% class(ms)) 
        ms <- NULL
    if (!is.null(ms)) 
        cmd <- get.meme.consensus(cmd, ...)
    if (grepl("-cons $none", cmd, fixed = T)) 
        cmd <- gsub("-cons $none", "", cmd, fixed = T)
    bg.list <- genome.info$bg.list[[seq.type]]
    bg.fname <- genome.info$bg.fname[seq.type]
    bgo <- bg.order[seq.type]
    if (TRUE && !force && !is.null(ms) && !is.null(ms$prev.run) && 
        length(ms$prev.run$rows) == length(rows) && all(ms$prev.run$rows %in% 
        rows) && all(rows %in% ms$prev.run$rows) && cmd == ms$prev.run$cmd && 
        (ms$prev.run$bg.order == bgo || sum(is.na(c(ms$prev.run$bg.order, 
            bgo) == 1))) && all(motif.upstream.scan[[seq.type]] == 
        ms$prev.run$m.u.scan)) {
        message("SKIPPING CLUSTER (UNCHANGED): ", k)
        ms$iter = iter
        ms$last.run = TRUE
        return(ms)
    }
    seqs <- get.sequences(rows, seq.type = seq.type, verbose = verbose, 
        ...)
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k, last.run = FALSE))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k, last.run = FALSE))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    all.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    if (is.null(all.seqs)) 
        all.seqs <- get.sequences("all", seq.type = seq.type, 
            distance = motif.upstream.scan[[seq.type]], filter = F, 
            ...)
    if (!is.na(bgo) && grepl("-bfile $bgFname", cmd, fixed = T)) {
        if (!recalc.bg && !file.exists(bg.fname)) {
            genome.info$bg.fname[seq.type] <- bg.fname <- my.tempfile("meme.tmp", 
                suf = ".bg")
            capture.output(genome.info$bg.list[[seq.type]] <- mkBgFile(all.seqs, 
                order = bgo, bgfname = bg.fname, input.list = bg.list, 
                use.rev.comp = grepl("-revcomp", meme.cmd[seq.type])))
        }
        if (recalc.bg || is.null(bg.list)) {
            tmp.seqs <- all.seqs[!names(all.seqs) %in% rows]
            bg.fname <- my.tempfile("meme.tmp", suf = ".bg")
            capture.output(bg.list <- mkBgFile(tmp.seqs, order = bg.order[seq.type], 
                verbose = F, bgfname = bg.fname, use.rev.comp = grepl("-revcomp", 
                  meme.cmd[seq.type])))
            rm(tmp.seqs)
        }
    }
    if (is.na(bgo) || is.null(bg.list)) {
        cmd <- gsub("-bfile $bgFname", "", cmd, fixed = T)
        bg.list <- NULL
    }
    if (exists("get.sequence.psps")) 
        psps <- get.sequence.psps(seqs)
    else psps <- NULL
    run.meme <- function(sgenes, seqs, cmd, seq.type, ...) {
        if (pal.opt == "non") {
            out <- runMeme(sgenes, seqs, cmd, seq.type = seq.type, 
                ...)
        }
        else if (pal.opt == "pal") {
            cmd <- paste(cmd, "-pal")
            out <- runMeme(sgenes, seqs, cmd, seq.type = seq.type, 
                ...)
        }
        else if (pal.opt == "both") {
            cmd2 <- paste(cmd, "-pal")
            out <- list(non = runMeme(sgenes, seqs, cmd, seq.type = seq.type, 
                ...), pal = runMeme(sgenes, seqs, cmd2, seq.type = seq.type, 
                ...))
        }
        out
    }
    if (verbose) {
        meme.out <- run.meme(names(seqs), seqs, cmd, verbose = verbose, 
            bg.list = bg.list, bgfname = bg.fname, psps = psps, 
            seq.type = seq.type, ...)
    }
    else {
        capture.output(meme.out <- try(run.meme(names(seqs), 
            seqs, cmd, verbose = verbose, bg.list = bg.list, 
            bgfname = bg.fname, psps = psps, seq.type = seq.type, 
            ...)))
    }
    if (pal.opt == "both") {
        meme.out2 <- lapply(meme.out, getMemeMotifInfo)
        e.vals <- sapply(lapply(meme.out2, sapply, "[[", "e.value"), 
            min, na.rm = T)
        e.vals[1] <- sqrt(e.vals[1])
        if (verbose) 
            cat(k, "Using pal/non-pal motif:", e.vals, names(which.min(e.vals)), 
                "\n")
        meme.out <- meme.out[[which.min(e.vals)]]
        meme.out2 <- meme.out2[[which.min(e.vals)]]
        attr(meme.out2, "is.pal") <- names(which.min(e.vals)) == 
            "pal"
    }
    else {
        meme.out2 <- getMemeMotifInfo(meme.out)
        attr(meme.out2, "meme.command.line") <- attr(meme.out, 
            "meme.command.line")
        attr(meme.out2, "is.pal") <- pal.opt == "pal"
    }
    if (length(meme.out2) <= 0) 
        return(list(k = k, last.run = FALSE))
    if (is.na(bgo) || is.null(bg.list) || !file.exists(bg.fname)) 
        bg.fname <- NULL
    if (verbose) 
        mast.out <- try(runMast(meme.out, mast.cmd[seq.type], 
            names(all.seqs), all.seqs, verbose = verbose, seq.type = seq.type, 
            bg.list = bg.list, bgfname = bg.fname, ...))
    else capture.output(mast.out <- try(runMast(meme.out, mast.cmd[seq.type], 
        names(all.seqs), all.seqs, verbose = verbose, seq.type = seq.type, 
        bg.list = bg.list, bgfname = bg.fname, ...)))
    pv.ev <- get.pv.ev.single(mast.out, rows)
    if (recalc.bg && !is.null(bg.fname) && file.exists(bg.fname) && 
        !"unlink" %in% names(list(...))) 
        unlink(bg.fname)
    prev.run <- list(rows = rows, cmd = cmd, bg.order = bgo, 
        m.u.scan = motif.upstream.scan[[seq.type]])
    invisible(list(k = k, last.run = FALSE, meme.out = meme.out2, 
        pv.ev = pv.ev, prev.run = prev.run))
}
mkBgFile <-
function (bgseqs = NULL, order = 0, bgfname = NULL, input.list = NULL, 
    use.rev.comp = T, verbose = T) 
{
    if (!is.null(input.list) && !is.null(bgfname)) {
        tmp <- unlist(input.list[2:length(input.list)])
        tmp2 <- sprintf("%.8f", tmp)
        names(tmp2) <- names(tmp)
        write.table(tmp2, row.names = names(tmp2), col.names = paste("#", 
            order, "th order Markov background model"), quote = F, 
            file = bgfname)
        return(input.list)
    }
    repl <- list(R = c("G", "A"), Y = c("T", "C"), K = c("G", 
        "T"), M = c("A", "C"), S = c("G", "C"), W = c("A", "T"), 
        N = c("G", "A", "T", "C"))
    bad.seqs <- grep("[^GATCX]", bgseqs, perl = T)
    if (length(bad.seqs) > 0) {
        if (verbose) 
            message(length(bad.seqs), " sequences with degenerate residues...fixing.")
        for (i in bad.seqs) {
            tmp <- strsplit(bgseqs[i], character(0))[[1]]
            inds <- grep("[^GATCX]", tmp, perl = T)
            for (ind in inds) tmp[ind] <- sample(repl[[tmp[ind]]], 
                1)
            bgseqs[i] <- paste(tmp, collapse = "")
        }
    }
    if (verbose) 
        cat("Calculating", order, "th order background Markov model from", 
            length(bgseqs), "sequences\n")
    if (use.rev.comp && verbose) 
        cat("Using reverse-complement too.\n")
    if (use.rev.comp) 
        bgseqs <- c(bgseqs, rev.comp(bgseqs))
    bgseqs <- bgseqs[!get.dup.seqs(bgseqs)]
    mc <- get.parallel(order + 1)
    apply.func <- lapply
    tmp <- mc$apply(0:order, function(ord, mc.cores) {
        out <- list()
        if (verbose) 
            cat("Calculating", ord, "th order part of background Markov model from", 
                length(bgseqs), "sequences\n")
        if (ord == 0) {
            all.substrings <- unlist(strsplit(bgseqs, character(0)), 
                use.names = F)
        }
        else {
            all.substrings <- sapply(1:(max(nchar(bgseqs)) - 
                ord), function(i) substr(bgseqs, i, i + ord))
            all.substrings <- as.vector(all.substrings)
        }
        all.substrings <- all.substrings[!is.na(all.substrings) & 
            all.substrings != "" & nchar(all.substrings) == ord + 
            1]
        counts <- table(as.factor(all.substrings))
        counts <- sort(counts)
        counts <- counts/length(all.substrings)
        counts <- counts[grep("N", names(counts), val = T, invert = T)]
        out <- as.list(counts)
        for (i in names(out)) {
            names(out[[i]]) <- NULL
            if (verbose && ord <= 3) 
                cat("FREQ:", i, "=", counts[i], "\n")
        }
        out
    }, mc.cores = min(order + 1, mc$par))
    out <- list()
    out$order <- order
    for (i in 1:length(tmp)) for (j in 1:length(tmp[[i]])) out[[names(tmp[[i]])[j]]] <- tmp[[i]][[j]]
    if (!is.null(bgfname) && !file.exists(bgfname)) {
        cat("Writing to file:", bgfname, "\n")
        tmp <- unlist(out)
        tmp <- tmp[2:length(tmp)]
        tmp2 <- sprintf("%.8f", tmp)
        names(tmp2) <- names(out)[2:length(out)]
        write.table(tmp2, row.names = names(tmp2), col.names = paste("#", 
            order, "th order Markov background model"), quote = F, 
            file = bgfname)
    }
    invisible(out)
}
mkTempMemeFiles <-
function (sgenes, seqs, fname = "meme.tmp.fst", bgseqs = NULL, 
    bgfname = NULL, filter.seqs = T, bg.list = NULL, force.overwrite = F, 
    seq.type = names(mot.weights)[1], seq.weights = NULL, psps = NULL, 
    ...) 
{
    if (!file.exists(fname) || file.info(fname)$size == 0 || 
        force.overwrite) {
        sgenes <- sgenes[!(is.na(seqs) | is.null(seqs) | seqs == 
            "")]
        seqs <- seqs[!(is.na(seqs) | is.null(seqs) | seqs == 
            "")]
        max.width <- as.integer(strsplit(meme.cmd[seq.type], 
            " ")[[1]][which(strsplit(meme.cmd[seq.type], " ")[[1]] == 
            "-maxw") + 1])
        if (filter.seqs) {
            sgenes <- sgenes[nchar(seqs) >= max.width]
            seqs <- seqs[nchar(seqs) >= max.width]
        }
        lengths <- sum(nchar(seqs)) + length(seqs) * 3
        if (!is.null(seq.weights)) {
            seq.weights <- seq.weights[sgenes]
            seq.weights[is.na(seq.weights)] <- 0
            cat(paste(">WEIGHTS", paste(seq.weights, collapse = " ")), 
                paste(">", sgenes, "\n", seqs, sep = ""), file = fname, 
                sep = "\n")
        }
        else {
            cat(paste(">", sgenes, "\n", seqs, sep = ""), file = fname, 
                sep = "\n")
        }
        if (!is.null(psps)) {
            psps <- psps[sgenes]
            for (i in names(psps)) {
                psps[[i]][is.na(psps[[i]]) | is.infinite(psps[[i]])] <- 0
                psps[[i]] <- psps[[i]]/(sum(psps[[i]], na.rm = T) + 
                  0.01)
            }
            psps <- sapply(psps[sgenes], function(i) paste(sprintf("%.5f", 
                i), collapse = " "))
            cat(paste(">", sgenes, "\n", psps, sep = ""), file = paste(fname, 
                ".psp", sep = ""), sep = "\n")
        }
    }
    if (force.overwrite || (!is.null(bgfname) && (!file.exists(bgfname) || 
        file.info(bgfname)$size <= 0))) {
        if (!is.null(bg.list)) 
            mkBgFile(input.list = bg.list, order = bg.list$order, 
                bgfname = bgfname)
        else if (!is.null(bgseqs)) 
            mkBgFile(bgseqs, order = 0, bgfname = bgfname)
    }
    length(seqs)
}
motif.all.clusters <-
function (ks = 1:k.clust, seq.type = names(mot.weights)[1], verbose = T, 
    debug = F, ...) 
{
    out.ms <- meme.scores[[seq.type]]
    mc <- get.parallel(length(ks), verbose = T, para.cores = get("parallel.cores.motif"))
    if (any(grepl("foreach", deparse(mc$apply))) && getDoParName() == 
        "doMC") 
        mc$apply <- function(list, FUN, ...) foreach(l = list, 
            .options.multicore = list(preschedule = F, set.seed = T)) %dopar% 
            {
                FUN(l, ...)
            }
    if (!debug) {
        out.ms <- mc$apply(ks, FUN = function(k) try(motif.one.cluster(k, 
            seq.type = seq.type, verbose = F, ...)))
    }
    else {
        message("DEBUG MODE: NOT PARALLELIZING!\n")
        out.ms <- lapply(ks, FUN = function(k) motif.one.cluster(k, 
            seq.type = seq.type, verbose = T, ...))
    }
    out.ms[[k.clust + 1]] <- ""
    for (k in ks) {
        if (length(out.ms) < k || is.null(out.ms[[k]]) || class(out.ms[[k]]) == 
            "try-error" || out.ms[[k]]$k != k || (!is.null(out.ms[[k]]$iter) && 
            out.ms[[k]]$iter != iter)) {
            out <- try(motif.one.cluster(k, seq.type = seq.type, 
                verbose = T, ...))
        }
        else {
            out <- out.ms[[k]]
        }
        if (class(out) == "try-error") 
            out <- try(motif.one.cluster(k, seq.type = seq.type, 
                verbose = T, ...))
        if (class(out) == "try-error" || is.null(out) || out$k != 
            k) {
            message("ERROR on cluster ", k)
            out <- list()
        }
        else if (verbose) {
            cat(iter, k, length(get.rows(k)), seq.type, "\t")
        }
        if (verbose) {
            if (is.null(out) || is.null(out$meme.out)) 
                cat("Inf \n")
            else {
                ind <- 1
                if (!is.null(out$pv.ev)) {
                  if ("p.value" %in% colnames(out$pv.ev[[1]])) 
                    mn <- mean(log10(out$pv.ev[[ind]][rownames(out$pv.ev[[ind]]) %in% 
                      get.rows(k), "p.value"]), na.rm = T)
                  else mn <- mean(log10(out$pv.ev[[ind]]$pvals), 
                    na.rm = T)
                }
                else {
                  mn <- "Inf"
                }
                cat(k, if (attr(out$meme.out, "is.pal")) 
                  "pal"
                else "non", sapply(out$meme.out[1:min(3, length(out$meme.out))], 
                  "[[", "e.value"), mn, "\t", pssm.to.string(out$meme.out[[1]]$pssm), 
                  "\n")
            }
        }
        out$iter <- iter
        out$k <- k
        out.ms[[k]] <- out
    }
    out.ms$all.pv <- make.pv.ev.matrix(out.ms)
    if (FALSE) {
        for (k in 1:k.clust) {
            m <- out.ms[[k]]
            if (!is.null(m) && !is.null(m$pv.ev)) 
                out.ms[[k]]$pv.ev[[1]] <- NULL
        }
    }
    attr(out.ms, "seq.type") <- seq.type
    invisible(out.ms)
}
motif.one.cluster <-
function (k, seq.type = names(mot.weights)[1], verbose = F, force = F, 
    ...) 
{
    st <- strsplit(seq.type, " ")[[1]]
    out <- meme.scores[[seq.type]][[k]]
    if (st[2] == "meme") {
        out <- meme.one.cluster(k, seq.type = seq.type, verbose, 
            force = force, ...)
    }
    else {
        if (is.numeric(k)) 
            rows <- get.rows(k)
        else rows <- k
        if (TRUE && !force) {
            if (!is.null(out) && !is.null(out$prev.run) && length(out$prev.run$rows) == 
                length(rows) && all(out$prev.run$rows %in% rows) && 
                all(rows %in% out$prev.run$rows) && all(motif.upstream.scan[[seq.type]] == 
                out$prev.run$m.u.scan)) {
                message("SKIPPING CLUSTER (UNCHANGED): ", k)
                out$iter = iter
                out$last.run = TRUE
                return(out)
            }
        }
        if (st[2] == "memepal") 
            out <- meme.one.cluster(k, seq.type = seq.type, verbose, 
                pal.opt = "pal", force = force, ...)
        else if (st[2] == "memeboth") 
            out <- meme.one.cluster(k, seq.type = seq.type, verbose, 
                pal.opt = "both", force = force, ...)
        else if (st[2] == "weeder") 
            out <- weeder.one.cluster(k, seq.type = seq.type, 
                verbose = verbose, n.motifs = 5, ...)
        else if (st[2] %in% c("spacer", "prism")) 
            out <- spacer.one.cluster(k, seq.type = seq.type, 
                verbose = verbose, ...)
        else if (st[2] == "cosmo") 
            out <- cosmo.one.cluster(k, seq.type = seq.type, 
                verbose = verbose, n.motifs = 1, ...)
        prev.run <- list(rows = rows, m.u.scan = motif.upstream.scan[[seq.type]])
        out$prev.run <- prev.run
    }
    invisible(out)
}
motif.similarities.custom <-
function (query = 1:(k.clust - 1), target = query, seq.type = names(mot.weights)[1], 
    e.value.cutoff = 100, resid.cutoff = 0.99, dist.meth = "ed", 
    min.overlap = 4, p.values = T, p.cutoff = 0.01, p.correct = T, 
    out.p.only = F, verbose = F, consensus = F, filter = T, add.to = NULL) 
{
    mc <- get.parallel(k.clust)
    mc$apply <- lapply
    query <- sort(query)
    target <- sort(target)
    outA <- do.call(c, mc$apply(query, function(i) {
        out <- list()
        if (verbose) 
            cat(i, "\n")
        if (length(get.rows(i)) <= 1) 
            return(out)
        if (clusterStack[[i]]$resid > resid.cutoff) 
            return(out)
        memeOut <- meme.scores[[seq.type]][[i]]$meme.out
        if (is.null(memeOut)) 
            return(out)
        for (ii in 1:length(memeOut)) {
            if (is.na(memeOut[[ii]]$e.value) || memeOut[[ii]]$e.value > 
                e.value.cutoff) 
                next
            pssm1 <- memeOut[[ii]]$pssm
            if (is.null(pssm1)) 
                next
            for (j in target[target >= i]) {
                if (length(get.rows(j)) <= 1) 
                  next
                if (clusterStack[[j]]$resid > resid.cutoff) 
                  next
                memeOut2 <- meme.scores[[seq.type]][[j]]$meme.out
                if (is.null(memeOut2)) 
                  next
                for (jj in 1:length(memeOut2)) {
                  if (is.na(memeOut2[[jj]]$e.value) || memeOut2[[jj]]$e.value > 
                    e.value.cutoff) 
                    next
                  pssm2 <- memeOut2[[jj]]$pssm
                  if (is.null(pssm2)) 
                    next
                  if (!is.null(add.to) && nrow(subset(add.to, 
                    biclust1 == i & motif1 == ii & biclust2 == 
                      j & motif2 == jj)) >= 1) 
                    next
                  tmp <- compare.pssms(pssm1, pssm2, rev = F, 
                    weight = F, score = dist.meth, min.ov = min.overlap)
                  out[[paste(i, ii, j, jj)]] <- cbind(i, ii, 
                    clusterStack[[i]]$resid, memeOut[[ii]]$e.value, 
                    memeOut[[ii]]$width, j, jj, clusterStack[[j]]$resid, 
                    memeOut[[jj]]$e.value, memeOut[[jj]]$width, 
                    tmp, nrow(tmp))
                  tmp <- compare.pssms(pssm1, pssm2, rev = T, 
                    weight = F, score = dist.meth, min.ov = min.overlap)
                  out[[paste(i, ii, j, -jj)]] <- cbind(i, ii, 
                    clusterStack[[i]]$resid, memeOut[[ii]]$e.value, 
                    memeOut[[ii]]$width, j, -jj, clusterStack[[j]]$resid, 
                    memeOut[[jj]]$e.value, memeOut[[ii]]$width, 
                    tmp, nrow(tmp))
                }
            }
        }
        out
    }))
    outA <- do.call(rbind, outA)
    rownames(outA) <- NULL
    if (is.null(outA)) 
        return(outA)
    colnames(outA) <- c("biclust1", "motif1", "resid1", "e.value1", 
        "width1", "biclust2", "motif2", "resid2", "e.value2", 
        "width2", "offset", "overlap", dist.meth, "n.tests")
    quantiles <- cdfs <- NULL
    out2 <- NULL
    ovs <- table(outA[, "overlap"])
    if (!is.na(p.cutoff)) {
        quantiles <- do.call(rbind, mc$apply(as.integer(names(ovs)), 
            function(ov) c(ov, quantile(outA[outA[, "overlap"] %in% 
                (ov + (if (ovs[as.character(ov)] < 1000) c(-1, 
                  0, 1) else 0)), dist.meth], 1 - p.cutoff))))
    }
    if (p.values) {
        cdfz <- mc$apply(as.integer(names(ovs)), function(ov) ecdf(outA[outA[, 
            "overlap"] %in% (ov + (if (ovs[as.character(ov)] < 
            1000) 
            c(-1, 0, 1)
        else 0)), dist.meth]))
        cdfs <- list()
        for (i in 1:length(ovs)) cdfs[[as.integer(names(ovs))[i]]] <- cdfz[[i]]
        rm(cdfz)
        qq <- numeric()
        qq[quantiles[, 1]] <- quantiles[, 2]
        print(dim(outA))
        if (filter) {
            out2 <- outA[outA[, dist.meth] > qq[outA[, "overlap"]], 
                ]
            rm(qq)
        }
        else out2 <- outA
        if (out.p.only) 
            rm(outA)
        print(dim(out2))
        pvs <- rep(NA, nrow(out2))
        tmp <- sapply(as.integer(names(ovs)), function(i) {
            pvs[which(out2[, "overlap"] == i)] <- cdfs[[i]](out2[which(out2[, 
                "overlap"] == i), dist.meth])
            pvs
        })
        if (out.p.only) 
            rm(cdfs)
        print(dim(out2))
        pvs <- 1 - apply(tmp, 1, sum, na.rm = T)
        rm(tmp)
        pvs <- pvs * out2[, "n.tests"]/2/2
        out2 <- cbind(out2, p.value = pvs)
        print(dim(out2))
        if (filter && !is.na(p.cutoff)) {
            out2 <- out2[pvs <= p.cutoff, ]
            pvs <- pvs[pvs <= p.cutoff]
        }
        out2 <- out2[order(pvs), ]
        print(dim(out2))
        rownames(out2) <- NULL
        colnames(out2)[1:14] <- c("biclust1", "motif1", "resid1", 
            "e.value1", "width1", "biclust2", "motif2", "resid2", 
            "e.value2", "width2", "offset", "overlap", dist.meth, 
            "n.tests")
    }
    if (consensus) {
        out2 <- as.data.frame(out2)
        for (ind in c("biclust1", "motif1", "width1", "biclust2", 
            "motif2", "width2", "offset", "overlap", "n.tests")) out2[[ind]] <- as.integer(out2[[ind]])
        out2 <- cbind(out2, as.data.frame(t(apply(out2, 1, function(i) c(pssm.to.string(meme.scores[[seq.type]][[i["biclust1"]]]$meme.out[[i["motif1"]]]$pssm), 
            pssm.to.string(meme.scores[[seq.type]][[i["biclust2"]]]$meme.out[[i["motif2"]]]$pssm))))))
        colnames(out2)[c(-1, 0) + ncol(out2)] <- c("consensus1", 
            "consensus2")
    }
    print(dim(out2))
    if (!out.p.only) 
        return(list(out = as.data.frame(outA), out.p = as.data.frame(out2), 
            quantiles = quantiles, cdfs = cdfs))
    else return(as.data.frame(out2))
}
motif.similarities.tomtom <-
function (query = 1:k.clust, target = 1:k.clust, query.mot = NA, 
    target.mot = NA, seq.type = names(mot.weights)[1], e.value.cutoff = 100, 
    resid.cutoff = 0.8, dist.meth = "ed", q.thresh = 0.5, min.overlap = 4, 
    q.pseudo = 0, t.pseudo = 0, min.gene.overlap = NA, desymmetrize = T, 
    unlink = T, files.only = F, verbose = T, ...) 
{
    meme.let <- c("A", "C", "G", "T")
    lines <- c("MEME version 3.0", "", "ALPHABET= ACGT", "", 
        "strands: + -", "", "Background letter frequencies (from dataset with add-one prior applied):")
    lines <- c(lines, paste(names(unlist(genome.info$bg.list[[seq.type]][meme.let])), 
        sprintf("%.3f", unlist(genome.info$bg.list[[seq.type]][meme.let])), 
        collapse = " "))
    query <- query[!is.na(query)]
    target <- target[!is.na(target)]
    if (is.na(query.mot)) 
        query.mot <- rep(NA, length(query))
    if (is.na(target.mot)) 
        target.mot <- rep(NA, length(target))
    cluster.motif.lines <- function(k, mot) {
        lines <- character()
        memeOut <- meme.scores[[seq.type]][[k]]
        if (is.null(memeOut) || memeOut == "") 
            return(lines)
        memeOut <- meme.scores[[seq.type]][[k]]$meme.out
        if (is.null(memeOut)) 
            return(lines)
        if (clusterStack[[k]]$resid > resid.cutoff) 
            return(lines)
        for (i in 1:length(memeOut)) {
            if (!is.na(mot) && i != mot) 
                next
            if (memeOut[[i]]$e.value > e.value.cutoff) 
                next
            pssm <- memeOut[[i]]$pssm
            lines <- c(lines, "", sprintf("MOTIF bic_%03d_%02d_%.3f_%.3e", 
                k, i, clusterStack[[k]]$resid, memeOut[[i]]$e.value), 
                sprintf("BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", 
                  k, i, clusterStack[[k]]$resid, memeOut[[i]]$e.value), 
                sprintf("letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", 
                  nrow(pssm), memeOut[[i]]$sites, memeOut[[i]]$e.value))
            lines <- c(lines, apply(pssm, 1, function(i) sprintf("%5.3f %5.3f %5.3f %5.3f", 
                i[1], i[2], i[3], i[4])))
        }
        lines
    }
    cmd <- "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"
    if (is.na(min.gene.overlap)) {
        lines.t <- c(lines, do.call(c, lapply(target, function(k) cluster.motif.lines(k, 
            target.mot[which(target == k)]))))
        if (verbose) 
            cat("TARGET MOTIFS:", length(grep("MOTIF", lines.t))/2, 
                range(target), "\n")
        tfile <- my.tempfile("tomtom_t_", )
        cat(lines.t, file = tfile, sep = "\n")
        cmd <- sprintf(cmd, progs.dir, q.thresh, dist.meth, min.overlap, 
            q.pseudo, t.pseudo, tfile)
    }
    if (!is.na(min.gene.overlap)) 
        c.rows <- lapply(target, get.rows)
    mc <- get.parallel(length(query))
    if (is.na(mc$par)) 
        mc$par <- 1
    if (mc$par > length(query)) 
        mc$par <- length(query)
    if (files.only == TRUE || (!is.logical(files.only) && !is.na(files.only) && 
        !is.null(files.only))) {
        mc$apply <- lapply
        mc$par <- 1
    }
    tout <- do.call(rbind, mc$apply(1:mc$par, function(par) {
        ks <- query[seq(par, length(query), by = mc$par)]
        lines.q <- c(lines, do.call(c, lapply(ks, function(k) cluster.motif.lines(k, 
            query.mot[which(query == k)]))))
        if (is.na(query.mot)) 
            n.query.motifs <- length(grep("MOTIF", lines.q))/2
        else n.query.motifs <- length(grep("MOTIF", lines.q))
        if (verbose) 
            cat("QUERY MOTIFS:", n.query.motifs, range(ks), "\n")
        if (n.query.motifs <= 0) 
            return(NULL)
        rows <- NULL
        if (!is.na(min.gene.overlap)) {
            rows <- get.rows(k)
            t.ok <- sapply(c.rows, function(r) sum(r %in% rows) >= 
                min.gene.overlap)
            t.ok[k] <- FALSE
            lines.t <- c(lines, do.call(c, lapply(target[t.ok], 
                function(k) cluster.motif.lines(k, target.mot[which(target == 
                  k)]))))
            if (verbose) 
                cat("TARGET MOTIFS:", length(grep("MOTIF", lines.t))/2, 
                  range(target), "\n")
            tfile <- my.tempfile("tomtom_t_", )
            cat(lines.t, file = tfile, sep = "\n")
            cmd <- sprintf(cmd, progs.dir, q.thresh, dist.meth, 
                min.overlap, q.pseudo, t.pseudo, tfile)
        }
        qfile <- paste(gsub("tomtom_t_", "tomtom_q_", tfile), 
            "_", min(ks, na.rm = T), "_", max(ks, na.rm = T), 
            sep = "")
        cat(lines.q, file = qfile, sep = "\n")
        cmd <- paste(cmd, "-query", qfile)
        if (files.only == TRUE || (!is.logical(files.only) && 
            !is.na(files.only) && !is.null(files.only))) {
            if (is.character(files.only)) 
                cat(paste(cmd, " > ", qfile, ".out\n", sep = ""), 
                  file = files.only, append = T)
            return(paste(cmd, " > ", qfile, ".out", sep = ""))
        }
        if (verbose) 
            print(cmd)
        tout <- system(cmd, intern = T)
        if (unlink) 
            unlink(qfile)
        tout <- do.call(rbind, strsplit(tout, "\t"))
        if (!is.null(tout)) {
            colnames(tout) <- tout[1, , drop = F]
            tout <- tout[-1, , drop = F]
            tout <- as.data.frame(tout[tout[, 1] != tout[, 2], 
                , drop = F])
        }
        print(dim(tout))
        tout
    }))
    if (files.only == TRUE || (!is.logical(files.only) && !is.na(files.only) && 
        !is.null(files.only))) 
        return(tout)
    cat("GOT", nrow(tout), "motif alignments.\n")
    q.id <- t(sapply(strsplit(as.character(tout[, 1]), "_"), 
        function(i) as.integer(i[2:3])))
    q.res.ev <- t(sapply(strsplit(as.character(tout[, 1]), "_"), 
        function(i) as.numeric(i[4:5])))
    t.id <- t(sapply(strsplit(as.character(tout[, 2]), "_"), 
        function(i) as.integer(i[2:3])))
    t.res.ev <- t(sapply(strsplit(as.character(tout[, 2]), "_"), 
        function(i) as.numeric(i[4:5])))
    tout2 <- data.frame(biclust1 = as.integer(q.id[, 1]), motif1 = as.integer(q.id[, 
        2]), resid1 = as.numeric(q.res.ev[, 1]), e.value1 = as.numeric(q.res.ev[, 
        2]), biclust2 = as.integer(t.id[, 1]), motif2 = as.integer(t.id[, 
        2]), resid2 = as.numeric(t.res.ev[, 1]), e.value2 = as.numeric(t.res.ev[, 
        2]), offset = as.integer(as.character(tout$`Optimal offset`)), 
        p.value = as.numeric(as.character(tout$`p-value`)), q.value = as.numeric(as.character(tout$`q-value`)), 
        overlap = as.integer(as.character(tout$Overlap)), consensus1 = as.factor(as.character(tout$`Query consensus`)), 
        consensus2 = as.factor(ifelse(tout$Orientation == "-", 
            rev.comp(as.character(tout$`Target consensus`)), 
            as.character(tout$`Target consensus`))), orientation = tout$Orientation)
    if (exists(cmd)) 
        attr(tout2, "tomtom.cmd") <- cmd
    tout2 <- tout2[order(tout2$q.value, tout2$p.value), ]
    rm(tout)
    if (desymmetrize) 
        tout2 <- desymmetrize.tomtom.results(tout2)
    tout2
}
my.tempfile <-
function (pattern = "file", tmpdir = tempdir(), suffix = "", 
    n.rnd.char = 20) 
{
    file.path(paste(tmpdir, "/", pattern, "_", paste(sample(c(LETTERS, 
        letters, 0:9, 0:9, 0:9, 0:9), n.rnd.char), collapse = ""), 
        suffix, sep = ""))
}
network <-
function (genes, conditions, q.val = T) 
{
    if (missing(genes)) 
        genes <- attr(e$ratios, "rnames")
    if (missing(conditions)) 
        conditions <- NULL
    all.rm <- get.biclusters(genes)
    names(all.rm) <- genes
    if (!is.null(conditions)) {
        all.cm <- get.biclusters(conditions)
        names(all.cm) <- conditions
    }
    require(multicore)
    apply.func <- mclapply
    out <- apply.func(genes, function(g) {
        print(g)
        tmp <- agglom(src = g, srcType = "gene", targetType = "gene", 
            path = "bicluster", q.val = q.val, cond.filter = conditions)
        done <- genes[1:which(genes == g)]
        tmp <- subset(tmp, p.value < 1 & count > 2 & !rownames(tmp) %in% 
            done)
        if (nrow(tmp) <= 0) 
            return(NULL)
        out <- data.frame(gene1 = g, int = "gene_gene", gene2 = rownames(tmp), 
            count = tmp[, 1], p.value = tmp[, 2], q.value = tmp[, 
                3])
        rownames(out) <- NULL
        print(dim(out))
        out
    })
    out <- do.call(rbind, out)
    out2 <- apply.func(1:length(tt.out2), function(mc) {
        mc <- paste("MOTC", mc, sep = "_")
        print(mc)
        tmp <- agglom(src = mc, srcType = "motif.cluster", targetType = "gene", 
            path = "motif", q.val = q.val, cond.filter = conditions)
        tmp <- subset(tmp, count > 2)
        if (nrow(tmp) <= 0) 
            return(NULL)
        out2 <- data.frame(gene1 = mc, int = "motc_gene", gene2 = rownames(tmp), 
            count = tmp[, 1], p.value = tmp[, 2], q.value = tmp[, 
                3])
        rownames(out2) <- NULL
        print(dim(out2))
        out2
    })
    out2 <- do.call(rbind, out2)
    noa <- data.frame()
    genes <- unique(c(as.character(subset(out, int == "gene_gene")$gene1), 
        as.character(subset(out, int == "gene_gene")$gene2)))
    noa <- rbind(noa, data.frame(name = genes, type = "gene"))
    motcs <- unique(as.character(subset(out2, int == "motc_gene")$gene1))
    noa <- rbind(noa, data.frame(name = motcs, type = "motif_cluster"))
    list(out = out, out2 = out2, noa = noa)
}
pareto.adjust.weights <-
function (kwin = 51, diter = 21, max.delta = 0.05) 
{
    out.scaling <- c(row = row.scaling[iter - 1], mot = mot.scaling[iter - 
        1], net = net.scaling[iter - 1])
    orig.scaling <- c(row = row.scaling[iter], mot = mot.scaling[iter], 
        net = net.scaling[iter])
    if (iter < diter/2) 
        return(out.scaling)
    diter <- min(diter, iter)
    kwin <- min(kwin, iter)
    avs <- NULL
    for (i in names(out.scaling)) {
        if (i == "row") 
            col <- "row.scores"
        else if (i == "mot") 
            col <- "mot.scores"
        else col <- paste(i, "scores", sep = ".")
        tmp <- stats[, col]
        tmp <- tmp[!is.na(tmp)]
        if (length(tmp) > 0) {
            av <- runmed(tmp, k = kwin)
            dy <- (av[length(av)] - av[max(1, length(av) - diter)])/diff(range(av, 
                na.rm = T))
            if (!is.na(dy)) 
                out.scaling[i] <- out.scaling[i] + min(0.1, dy)
            out.scaling[i] <- max(out.scaling[i], orig.scaling[i])
        }
    }
    out.scaling
}
pareto.adjust.weights.OLD <-
function (iter, delta.iter = 200, delta.factor = 1, n.avg = 50, 
    max.delta = 0.05) 
{
    if (iter == 1) 
        return(c(row = row.scaling[iter], mot = mot.scaling[iter], 
            net = net.scaling[iter]))
    out.scaling <- c(row = row.scaling[iter - 1], mot = mot.scaling[iter - 
        1], net = net.scaling[iter - 1])
    if (iter < delta.iter + n.avg + 10) 
        return(out.scaling)
    all.diffs <- numeric()
    for (i in c("row", "mot", "net")) {
        if (i == "row") 
            col <- "resid"
        else if (i == "mot") 
            col <- "p.clust"
        else col <- paste(i, "scores", sep = ".")
        tops <- stats[[col]][(iter - n.avg + 1):iter]
        bots <- stats[[col]][(iter - n.avg + 1):iter - delta.iter]
        all.diffs[i] <- mean(tops - bots, na.rm = T)/abs(diff(range(stats[[col]], 
            na.rm = T)))
    }
    for (i in which(all.diffs > 0 & !is.na(all.diffs) & !is.na(out.scaling) & 
        out.scaling > 0)) out.scaling[names(all.diffs)[i]] <- out.scaling[names(all.diffs)[i]] + 
        min(c(out.scaling[i] * max.delta, all.diffs[i] * delta.factor))
    out.scaling
}
parse.blast.out <-
function (blast.out) 
{
    if (substr(blast.out[1], 1, 12) == "# BLASTN 2.2") 
        blast.out <- blast.out[-(1:3)]
    out <- t(sapply(strsplit(blast.out, "\t"), cbind))
    out <- data.frame(`Query id` = out[, 1], `Subject id` = out[, 
        2], `% identity` = as.numeric(out[, 3]), `alignment length` = as.integer(out[, 
        4]), mismatches = as.integer(out[, 5]), `gap openings` = as.integer(out[, 
        6]), `q. start` = as.integer(out[, 7]), `q. end` = as.integer(out[, 
        8]), `s. start` = as.integer(out[, 9]), `s. end` = as.integer(out[, 
        10]), `e-value` = as.numeric(out[, 11]), `bit score` = as.numeric(out[, 
        12]))
    out
}
plot.cluster.scores <-
function (scores) 
{
    row.memb <- attr(ratios, "rnames") %in% get.rows(scores$k)
    names(row.memb) <- attr(ratios, "rnames")
    par(mfrow = c(2, 2))
    for (i in colnames(scores$r)[1:4]) {
        h <- hist(scores$r[, i], breaks = 200, main = i)
        hist(rep(scores$r[row.memb, i], 10), breaks = h$breaks, 
            add = T, col = "red", border = "red")
        p <- scores$r.d[, i]
        lines(sort(scores$r[, i]), p[order(scores$r[, i])]/max(p) * 
            max(h$counts)/2, col = "green")
    }
}
plotClust <-
function (k, cluster = NULL, w.motifs = T, all.conds = T, title = NULL, 
    o.genes = NULL, dont.plot = F, network = "all", short.names = organism == 
        "sce", seq.type = names(mot.weights), ...) 
{
    if (!dont.plot) 
        opar <- par(no.readonly = T)
    if (!is.null(cluster)) {
        if (!dont.plot) 
            plotCluster.motif(cluster, seqs = cluster$seqs, p.val.shade.cutoff = 1, 
                o.genes = o.genes, no.plotCluster = all.conds, 
                ...)
        return(invisible(cluster))
    }
    c <- get.clust(k, varNorm = F)
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
                if ("gene" %in% colnames(c[[st]]$motif.out$pv.ev[[1]])) 
                  c[[st]]$motif.out$pv.ev[[2]] <- c[[st]]$motif.out$pv.ev[[1]]
                if (!is.null(meme.scores[[st]]$all.pv)) {
                  tmp <- cbind(p.value = meme.scores[[st]]$all.pv[, 
                    k], e.value = if ("all.ev" %in% names(meme.scores[[st]])) 
                    meme.scores[[st]]$all.ev[, k]
                  else NA)
                }
                else {
                  pv.ev <- meme.scores[[st]][[k]]$pv.ev[[1]]
                  if (ncol(pv.ev) <= 2) 
                    pv.ev <- meme.scores[[st]][[k]]$pv.ev[[2]]
                  tmp <- as.matrix(pv.ev[, 2:ncol(pv.ev)])
                  rownames(tmp) <- pv.ev[, 1]
                  colnames(tmp) <- c("p.value", "posns", "mots")
                }
                c[[st]]$motif.out$pv.ev[[1]] <- tmp
                c[[st]]$motif.out$p.values <- log10(c[[st]]$motif.out$pv.ev[[1]][, 
                  "p.value"])
                names(c[[st]]$motif.out$p.values) <- rownames(c[[st]]$motif.out$pv.ev[[1]])
            }
        }
    }
    if (!is.na(mot.iters[1]) && !no.genome.info) {
        c$seqs <- get.sequences(rows, distance = motif.upstream.scan[[seq.type[1]]], 
            seq.type = seq.type[1], filter = T, uniq = F)
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
            if (!i %in% names(networks)) 
                next
            tmp.net <- networks[[i]][networks[[i]]$protein1 %in% 
                rows & networks[[i]]$protein2 %in% rows, ]
            tmp.net <- cbind(tmp.net, net = rep(i, nrow(tmp.net)))
            c$network <- if (!is.null(c$network)) 
                rbind(c$network, tmp.net)
            else tmp.net
        }
    }
    c$gene.coords <- get.long.names(rows, short = short.names)
    if ("cog.code" %in% names(genome.info)) 
        c$cog.code <- genome.info$cog.code[rows]
    if (!is.null(title)) 
        c$name <- title
    if (!dont.plot) {
        plotCluster.motif(c, seqs = c$seqs, p.val.shade.cutoff = 1, 
            o.genes = o.genes, no.plotCluster = all.conds, ...)
        if (!"layout" %in% names(list(...))) 
            par(opar)
    }
    invisible(c)
}
plotCluster <-
function (cluster, imag = F, cond.labels = F, o.genes = NULL, 
    col.func = if (imag) topo.colors else rainbow, rats.names = names(ratios), 
    main = NULL, range.r = NULL, no.par = F, sort = F, box.plot = F, 
    ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    k <- cluster$k
    if (is.null(main)) 
        main <- paste(sprintf("Cluster: %04d %s; resid: %s; r/c: %d/%d", 
            k, organism, paste(sprintf("%.2f", cluster$resid[rats.names]), 
                collapse = " "), length(cluster$rows), length(cluster$cols)))
    rats <- get.cluster.matrix(unique(c(cluster$rows, o.genes)), 
        cluster$cols, matrices = rats.names)
    cols.b <- colnames(rats)[colnames(rats) %in% cluster$cols]
    if (sort) {
        o1 <- order(apply(rats[cluster$rows, cols.b, drop = F], 
            2, mean, na.rm = T))
        cols.b <- cols.b[o1]
        rats <- rats[, cols.b, drop = F]
    }
    if (all(is.na(rats))) {
        plot(0, 0, typ = "n", min = main, ...)
        return()
    }
    if (is.vector(rats)) {
        rats <- t(rats)
        rownames(rats) <- cluster$rows
    }
    if (imag) {
        grey.image <- function(mat, n.gray = 32, x = 1:nrow(mat), 
            y = 1:ncol(mat), col = gray((0:n.gray)/n.gray), ...) image(x, 
            y, mat, col = col, ...)
        grey.image(t(rats), col = col.func(256))
        return()
    }
    if (is.null(range.r)) 
        range.r <- range(rats[rats != min(rats, na.rm = T) & 
            rats != max(rats, na.rm = T)], na.rm = T)
    if (cond.labels && cluster$ncols < 100) 
        range.r[1] <- range.r[1] * 1.5
    if (!no.par) 
        par(mar = rep(2, 4), mgp = c(3, 1, 0) * 0.5)
    plot(1:length(cols.b), ylim = range.r, xlab = NA, ylab = NA, 
        main = main, typ = "n", xaxs = "i", ...)
    if (length(rats.names) > 1) {
        ind <- 0.5
        rts <- NULL
        for (i in 1:length(rats.names)) {
            col <- sapply(col2rgb(i + 1)/255 + 0.9, function(cc) min(cc, 
                1))
            col <- rgb(col[1], col[2], col[3])
            rect(ind, range.r[1] + 0.05, ind + sum(colnames(rats) %in% 
                colnames(ratios[[rats.names[i]]])), range.r[2] - 
                0.05, col = col, dens = NA)
            ind <- ind + sum(colnames(rats) %in% colnames(ratios[[rats.names[i]]]))
            rts <- cbind(rts, rats[, cols.b[cols.b %in% colnames(ratios[[rats.names[i]]]), 
                drop = F]])
        }
        rats <- rts
        rm(rts)
        cols.b <- colnames(rats)
    }
    if (exists("col.rug")) {
        if (is.integer(col.rug)) 
            colmap <- col.func(max(col.rug))[col.rug[cols.b]]
        else colmap <- col.rug[cols.b]
    }
    else if (all(deparse(col.func) == deparse(rainbow))) {
        colmap <- col.func(length(cols.b))
    }
    else {
        colmap <- col.func(cols.b)
    }
    if (box.plot) {
        colMeans <- apply(rats[cluster$rows, , drop = F], 2, 
            mean, na.rm = T)
        colSd <- apply(rats[cluster$rows, , drop = F], 2, sd, 
            na.rm = T)
        matlines(1:length(cols.b), cbind(colMeans - 2 * colSd, 
            colMeans + 2 * colSd), lty = 1, col = "lightgrey")
        boxplot(as.data.frame(rats[cluster$rows, , drop = F]), 
            ylim = range.r, names = NA, main = main, col = colmap, 
            outline = FALSE, border = FALSE, add = T, xaxs = "i", 
            xaxt = "n", ...)
        if (sort) 
            lines(1:length(cols.b), colMeans, lty = 1, lwd = 1, 
                col = "red")
    }
    else {
        cmap <- col.func(cluster$nrows)
        matlines(1:length(cols.b), t(rats[cluster$rows, , drop = F]), 
            ylim = range.r, xlab = NA, ylab = NA, main = main, 
            col = cmap, lty = 1, ...)
        if (exists("col.rug")) 
            for (i in unique(col.rug)) rug(which(cols.b %in% 
                names(which(col.rug == i))), col = colmap[which(col.rug == 
                i)[1]])
    }
    if (cond.labels) {
        tmp.y <- rep(range.r[1] * 0.85, cluster$ncols)
        cols <- if (box.plot) 
            colmap
        else "black"
        text(1:cluster$ncols, tmp.y, cols.b, srt = 90, col = cols, 
            ...)
    }
    if (names(dev.cur()) == "devSVG") {
        par(family = "Arial")
        for (c in 1:length(cols.b)) {
            setSVGShapeToolTip(cols.b[c])
            rect(c, range.r[1], c + 1, range.r[2], col = NA, 
                border = NA)
        }
    }
    if (!is.null(o.genes)) {
        matlines(1:length(cols.b), t(rats[o.genes, , drop = F]), 
            lty = 1, lwd = 3, col = 2:6)
        legend("bottomright", legend = o.genes, lty = 1, lwd = 3, 
            col = 2:6, cex = 0.7, bty = "n")
    }
}
plotCluster.all.conds <-
function (cluster, imag = F, cond.labels = F, o.genes = NULL, 
    rats.names = names(ratios), range.r = NULL, sort = F, box.plot = F, 
    col.func = if (imag) topo.colors else rainbow, ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    k <- cluster$k
    main <- paste(sprintf("Cluster: %04d %s; resid: %s; r/c: %d/%d", 
        k, organism, paste(sprintf("%.2f", cluster$resid[rats.names]), 
            collapse = " "), length(cluster$rows), length(cluster$cols)))
    rats <- get.cluster.matrix(unique(c(cluster$rows, o.genes)), 
        NULL, matrices = rats.names)
    cols.b <- c(colnames(rats)[colnames(rats) %in% cluster$cols], 
        colnames(rats)[!colnames(rats) %in% cluster$cols])
    if (sort) {
        inClust <- colnames(rats)[colnames(rats) %in% cluster$cols]
        o1 <- order(apply(rats[cluster$rows, inClust, drop = F], 
            2, mean, na.rm = T))
        outClust <- colnames(rats)[!colnames(rats) %in% cluster$cols]
        o2 <- order(apply(rats[cluster$rows, outClust, drop = F], 
            2, mean, na.rm = T))
        cols.b <- c(inClust[o1], outClust[o2])
    }
    len.b <- length(cols.b)
    rats <- rats[, cols.b, drop = F]
    par(mar = rep(2, 4), mgp = c(3, 1, 0) * 0.5)
    if (all(is.na(rats))) {
        plot(0, 0, typ = "n", main = main, ...)
        return()
    }
    if (is.vector(rats)) {
        rats <- t(rats)
        rownames(rats) <- cluster$rows
    }
    if (imag) {
        grey.image(t(rats), col = col.func(256))
        lines(rep(cluster$ncols + 0.5, 2), c(-999, 9999), col = 2, 
            lwd = 3, lty = 2)
        return()
    }
    if (is.null(range.r)) 
        range.r <- range(rats[rats != min(rats, na.rm = T) & 
            rats != max(rats, na.rm = T)], na.rm = T)
    if (cond.labels && len.b < 100) 
        range.r[1] <- range.r[1] * 1.5
    plot(1:len.b, xlim = c(0.95, len.b + 0.05), ylim = range.r, 
        xlab = NA, ylab = NA, main = main, typ = "n", xaxs = "i", 
        ...)
    if (length(ratios) > 1) {
        ind <- 0.5
        rts.in <- rts.out <- NULL
        for (in.out in 1:2) {
            if (in.out == 1) 
                cols <- cols.b[cols.b %in% cluster$cols]
            else if (in.out == 2) 
                cols <- cols.b[!cols.b %in% cluster$cols]
            for (i in 1:length(ratios)) {
                col <- sapply(col2rgb(i + 1)/255 + 0.9, function(cc) min(cc, 
                  1))
                col <- rgb(col[1], col[2], col[3])
                rect(ind, range.r[1] + 0.01, ind + sum(cols %in% 
                  colnames(ratios[[i]])), range.r[2] - 0.01, 
                  col = col, dens = NA)
                ind <- ind + sum(cols %in% colnames(ratios[[i]]))
                if (in.out == 1) 
                  rts.in <- cbind(rts.in, rats[, cols[cols %in% 
                    colnames(ratios[[i]])], drop = F])
                else if (in.out == 2) 
                  rts.out <- cbind(rts.out, rats[, cols[cols %in% 
                    colnames(ratios[[i]])], drop = F])
            }
        }
        rats <- cbind(rts.in, rts.out)
        cols.b <- colnames(rats)
        rm(rts.in, rts.out)
    }
    if (exists("col.rug")) {
        if (is.integer(col.rug)) 
            colmap <- col.func(max(col.rug))[col.rug[cols.b]]
        else colmap <- col.rug[cols.b]
    }
    else if (all(deparse(col.func) == deparse(rainbow))) {
        colmap <- col.func(length(cols.b))
    }
    else {
        colmap <- col.func(cols.b)
    }
    if (box.plot) {
        colMeans <- apply(rats[cluster$rows, , drop = F], 2, 
            mean, na.rm = T)
        colSd <- apply(rats[cluster$rows, , drop = F], 2, sd, 
            na.rm = T)
        matlines(1:length(cols.b), cbind(colMeans - 2 * colSd, 
            colMeans + 2 * colSd), lty = 1, col = "lightgrey")
        boxplot(as.data.frame(rats[cluster$rows, , drop = F]), 
            ylim = range.r, names = NA, main = main, col = colmap, 
            outline = FALSE, border = FALSE, add = T, xaxs = "i", 
            xaxt = "n", ...)
        if (sort) 
            lines(1:length(cols.b), colMeans, lty = 1, lwd = 1, 
                col = "red")
    }
    else {
        cmap <- col.func(cluster$nrows)
        matlines(1:len.b, t(rats[cluster$rows, , drop = F]), 
            ylim = range.r, col = cmap, main = main, xlab = NA, 
            ylab = NA, lty = 1, ...)
        if (exists("col.rug")) 
            for (i in unique(col.rug)) rug(which(cols.b %in% 
                names(which(col.rug == i))), col = colmap[which(col.rug == 
                i)[1]])
    }
    cols.in <- colnames(rats)[colnames(rats) %in% cluster$cols]
    lines(rep(length(cols.in) + 0.5, 2), range.r, col = 2, lwd = 3, 
        lty = 2)
    if (!is.null(o.genes)) {
        matlines(1:len.b, t(rats[o.genes, , drop = F]), lty = 1, 
            lwd = 3, col = 2:6)
        legend("bottomright", legend = o.genes, lty = 1, lwd = 3, 
            col = 2:6, cex = 0.7, bty = "n")
    }
    if (cond.labels) {
        tmp.y <- rep(range.r[1] * 0.85, len.b)
        cols <- if (box.plot) 
            colmap
        else "black"
        text(1:len.b, tmp.y, cols.b, srt = 90, col = cols, ...)
    }
    if (names(dev.cur()) == "devSVG") {
        par(family = "Arial")
        for (c in 1:length(cols.b)) {
            setSVGShapeToolTip(cols.b[c])
            rect(c, range.r[1], c + 1, range.r[2], col = NA, 
                border = NA)
        }
    }
}
plotCluster.motif <-
function (cluster, seqs = cluster$seqs, layout = NULL, colors = NULL, 
    motif.e.cutoff = Inf, no.plotCluster = T, addl.text = NULL, 
    ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    if (names(dev.cur()) == "devSVG") 
        par(family = "Arial")
    if (any(!cluster$rows %in% attr(ratios, "rnames"))) {
        cluster$rows <- cluster$rows[cluster$rows %in% attr(ratios, 
            "rnames")]
        cluster$nrows <- length(cluster$rows)
        warning(cluster$k, ": Some cluster rows are not in the ratios. Will plot without these rows.\n")
    }
    if (any(!cluster$cols %in% attr(ratios, "cnames"))) {
        cluster$cols <- cluster$cols[cluster$cols %in% attr(ratios, 
            "cnames")]
        cluster$ncols <- length(cluster$cols)
        warning(cluster$k, ": Some cluster cols are not in the ratios. Will plot without these cols.\n")
    }
    seq.types <- cluster$seq.type
    if (length(seq.types) == 1) {
        n.pssm.plot <- 3
    }
    else {
        n.pssm.plot <- 6
    }
    if (is.null(layout)) {
        if (length(seq.types) == 1) {
            layout <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 8, 2, 
                2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 
                8, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 
                4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5, 5, 
                6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7), ncol = 17, 
                byrow = T)
        }
        else {
            layout <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 11, 2, 
                2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 
                11, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 
                1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 
                1, 1, 1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 
                3, 3, 4, 4, 4, 4, 10, 10, 10, 10, 10, 10, 10, 
                10, 10, 5, 5, 5, 5, 6, 6, 6, 6, 10, 10, 10, 10, 
                10, 10, 10, 10, 10, 7, 7, 7, 7, 9, 9, 9, 9, 10, 
                10, 10, 10, 10, 10, 10, 10, 10, 8, 8, 8, 8, 9, 
                9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10), 
                ncol = 17, byrow = T)
        }
        if (no.plotCluster) {
            layout[layout == 1 | layout == max(layout)] <- 2
            layout <- layout - 1
            layout[, 1][layout[, 1] == 1] <- max(layout) + 1
        }
    }
    layout(layout)
    k <- cluster$k
    if (!is.null(ratios)) {
        args <- list(...)
        args <- args[names(args) != "p.val.shade.cutoff"]
        args$cluster <- cluster
        do.call(plotCluster.all.conds, args)
        if (!no.plotCluster) 
            do.call(plotCluster, args)
    }
    rows <- cluster$rows
    if (is.null(colors) || !is.null(cluster$cog.code)) {
        tmp.lett <- 1:26
        names(tmp.lett) <- LETTERS
        if (!is.null(cluster$cog.code)) 
            coo <- cluster$cog.code[rows]
        else coo <- 1:length(rows)
        tmp <- unique(tmp.lett[coo])
        names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
        cols <- rainbow(length(tmp))
        names(cols) <- names(tmp)
        cols <- cols[names(tmp.lett[coo])]
        cols[is.na(names(cols))] <- "darkgrey"
        names(cols) <- rows
    }
    else {
        cols <- rainbow(length(rows))
        names(cols) <- rows
    }
    colors <- cols
    n.plotted <- 1
    for (seq.type in seq.types) {
        if (n.plotted > n.pssm.plot) 
            break
        if (is.null(cluster[[seq.type]]$e.val) || all(is.na(cluster[[seq.type]]$e.val)) || 
            is.null(cluster[[seq.type]]$motif.out) || is.null(cluster[[seq.type]]$motif.out$pssms)) 
            next
        pssm <- cluster[[seq.type]]$motif.out$pssms
        for (ppp in 1:min(floor(n.pssm.plot/length(seq.types)), 
            length(pssm))) {
            if (n.plotted > n.pssm.plot) 
                break
            if (cluster[[seq.type]]$motif.out$e.values[ppp] > 
                motif.e.cutoff) 
                next
            viewPssm(pssm[[ppp]], mot.ind = ppp, main.title = sprintf("%s PSSM #%d; e=%.3g", 
                seq.type, ppp, cluster[[seq.type]]$motif.out$e.values[ppp]), 
                cex.main = 0.9)
            n.plotted <- n.plotted + 1
        }
    }
    while (n.plotted <= n.pssm.plot) {
        plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
        n.plotted <- n.plotted + 1
    }
    suppressWarnings(cluster <- plotCluster.network(cluster, 
        ...))
    if (is.null(seqs)) {
        seqs <- rep("", length(cluster$rows))
        names(seqs) <- cluster$rows
    }
    if (!is.null(seq.type) && !is.null(seqs) && length(seqs) > 
        0) 
        plotClusterMotifPositions(cluster, seqs, colors = colors, 
            ...)
    else plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "")
    try({
        par(mar = rep(0.5, 4), mgp = c(2, 1, 0) * 0.5)
        plot(c(0.5, 2.5), c(-1, 1), type = "n", tck = 0.01, cex.lab = 0.2, 
            cex.sub = 0.2, cex.axis = 0.2, axes = F)
        if (names(dev.cur()) == "devSVG") {
            par(family = "Arial")
            setSVGShapeToolTip(title = paste("Cluster:", sprintf("%04d", 
                k), organism, cmonkey.version), desc1 = sprintf("resid = %s; genes = %d; conds = %d", 
                paste(sprintf("%.2f", cluster$resid), collapse = " "), 
                length(cluster$rows), length(cluster$cols)))
            setSVGShapeURL(paste("http://www.genome.ad.jp/dbget-bin/www_bget?", 
                paste(organism, ":", cluster$rows, sep = "", 
                  collapse = "+"), sep = ""))
            rect(0.5, -1, 3.25, +1, col = "lightgreen", border = NA)
        }
    })
    if (!is.null(cluster$name)) 
        text(0.65, 0, cluster$name, srt = 90, xpd = NA, cex = 1)
    else if (!is.null(addl.text)) 
        text(0.65, 0, addl.text, srt = 90, xpd = NA, cex = 1)
    text(1.5, 0, sprintf("%s iter=%d", date.run, iter), srt = 90, 
        xpd = NA, cex = 1)
    text(2.35, 0, paste("cMonkey Version", cmonkey.version, organism), 
        srt = 90, xpd = NA, cex = 1)
    invisible(cluster)
}
plotCluster.network <-
function (cluster, network = "all", o.genes = NULL, colors = NULL, 
    cex = 0.7, no.legend = F, ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    require(igraph)
    rows <- cluster$rows
    if (is.null(cluster$network)) {
        if (network == "all") 
            network <- names(networks)
        for (i in network) {
            if (!i %in% names(networks)) 
                next
            tmp.net <- networks[[i]][networks[[i]]$protein1 %in% 
                rows & networks[[i]]$protein2 %in% rows, ]
            tmp.net <- cbind(tmp.net, net = rep(i, nrow(tmp.net)))
            cluster$network <- if (!is.null(cluster$network)) 
                rbind(cluster$network, tmp.net)
            else tmp.net
        }
    }
    network <- cluster$network
    nrows <- rows
    if (!is.null(o.genes)) 
        nrows <- c(nrows, o.genes)
    if (is.null(cluster$cog.code) && "cog.code" %in% names(genome.info)) 
        cluster$cog.code <- genome.info$cog.code[rows]
    if (is.null(cluster$colors)) {
        if (is.null(colors) || "cog.code" %in% names(genome.info)) {
            tmp.lett <- 1:26
            names(tmp.lett) <- LETTERS
            if (!is.null(cluster$cog.code)) 
                coo <- cluster$cog.code[rows]
            else coo <- 1:length(rows)
            tmp <- unique(tmp.lett[coo])
            names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
            cols <- rainbow(length(tmp))
            names(cols) <- names(tmp)
            cols <- cols[names(tmp.lett[coo])]
            cols[is.na(names(cols))] <- "darkgrey"
            names(cols) <- rows
            cluster$colors <- cols
        }
        else {
            cols <- rainbow(length(rows))
            names(cols) <- rows
            cluster$colors <- cols
        }
    }
    colors <- cluster$colors
    if (is.null(network) || nrow(network) <= 0) 
        network <- data.frame(protein1 = nrows, protein2 = nrows, 
            combined_score = jitter(rep(1/50, length(nrows))), 
            net = rep("none", length(nrows)))
    not.in <- nrows[!nrows %in% network$protein1 & !nrows %in% 
        network$protein2]
    for (i in not.in) network <- rbind(network, data.frame(protein1 = i, 
        protein2 = i, combined_score = 0, net = "none"))
    gr <- graph.edgelist(as.matrix(network[, 1:2]), directed = F)
    net.wts <- as.numeric(network$combined_score)
    names(net.wts) <- as.character(network$net)
    for (n in unique(names(net.wts))) {
        if (n == "none") 
            next
        net.wts[names(net.wts) == n] <- net.wts[names(net.wts) == 
            n]/max(net.wts[names(net.wts) == n], na.rm = T)
    }
    gr.layout <- layout.fruchterman.reingold(gr, niter = 3000, 
        weights = net.wts/5)
    gr.layout <- layout.norm(gr.layout, -1, 1, -1, 1)
    edge.colors <- character()
    curves <- rep(0, nrow(network))
    nets <- unique(as.character(network$net))
    if ("none" %in% nets) {
        nets <- unique(c("none", nets))
        inds <- c(1, 1:(length(nets) - 1))
        inds <- inds[inds != 0]
        inds <- inds[1:length(nets)]
        net.colors <- t(col2rgb(inds, T))/255
        net.colors[1, 4] <- 0
    }
    else {
        net.colors <- t(col2rgb(1:length(nets), T))/255
    }
    rownames(net.colors) <- nets
    for (i in 1:nrow(network)) {
        net <- as.character(network$net)[i]
        shade <- net.wts[i]
        nodes <- c(as.character(network$protein1)[i], as.character(network$protein2)[i])
        sub.net <- subset(network, protein1 %in% nodes & protein2 %in% 
            nodes)
        sub.nets <- unique(as.character(sub.net$net))
        curve.it <- max(0, nrow(sub.net) - 2)/2 * 0.33
        col <- net.colors[net, ]
        col2 <- col
        col2[col2 == 0] <- 1 - shade
        edge.colors[i] <- if (names(dev.cur()) != "X11") 
            rgb(col[1], col[2], col[3], shade)
        else rgb(col2[1], col2[2], col2[3])
        curves[i] <- curve.it * floor(which(sub.nets == net)/2) * 
            (if (which(sub.nets == net)%%2 == 0) 
                -1
            else 1)
    }
    if (all(curves == 0)) 
        curves <- FALSE
    if (!no.legend) {
        labels <- try(get.long.names(get.vertex.attribute(gr, 
            "name"), short = T), silent = T)
        if (class(labels) == "try-error") 
            labels <- get.vertex.attribute(gr, "name")
        labels[is.na(labels) | labels == ""] <- get.vertex.attribute(gr, 
            "name")[is.na(labels) | labels == ""]
    }
    else labels <- NA
    plot(gr, layout = gr.layout, margin = 0, rescale = F, edge.curved = curves, 
        vertex.color = colors[get.vertex.attribute(gr, "name")], 
        vertex.frame.color = colors[get.vertex.attribute(gr, 
            "name")], vertex.label.cex = cex, vertex.size = 7, 
        vertex.label = labels, vertex.label.family = if (names(dev.cur()) != 
            "pdf") 
            "Arial"
        else "sans", edge.color = edge.colors, edge.width = round(net.wts) + 
            1)
    if (!no.legend && length(nets[nets != "none"]) > 0) 
        legend("bottomright", legend = nets[nets != "none"], 
            col = 1:length(nets[nets != "none"]), lty = 1, lwd = 2, 
            bty = "n", cex = 0.5)
    if (names(dev.cur()) == "devSVG") {
        names <- cluster$gene.coords
        for (i in 1:nrow(gr.layout)) {
            gene <- get.vertex.attribute(gr, "name")[i]
            setSVGShapeToolTip(title = gene, desc1 = ifelse(is.na(names[gene]), 
                "", names[gene]))
            setSVGShapeURL(paste("http://www.genome.ad.jp/dbget-bin/www_bget?", 
                organism, ":", gene, sep = ""))
            points(gr.layout[i, 1], gr.layout[i, 2], col = "#FF000001", 
                cex = 10/3)
        }
    }
    cluster
}
plotClusterMotifPositions <-
function (cluster, seqs = cluster$seqs, long.names = T, shade = T, 
    p.val.shade.cutoff = 999, colors = NULL, sort.by = "p.value", 
    o.genes = NULL, no.key = F, short.names = organism == "sce", 
    seq.type = cluster$seq.type[1], ...) 
{
    if (length(cluster$rows) <= 0) {
        warning("Trying to plot a cluster with no rows!")
        return()
    }
    k <- cluster$k
    rows <- cluster$rows
    if (!is.null(o.genes)) 
        rows <- c(rows, o.genes)
    motif.out <- NULL
    if (!is.null(seq.type) && seq.type %in% names(cluster)) 
        motif.out <- cluster[[seq.type]]$motif.out
    is.dup.seq <- get.dup.seqs(cluster$seqs)
    p.clust <- cluster$p.clust
    e.clust <- cluster$e.val
    motif.info <- NULL
    if ((!all(is.na(p.clust)) || !all(is.na(e.clust))) && !is.null(motif.out) && 
        !is.null(motif.out$pv.ev)) 
        motif.info <- subset(motif.out$pv.ev[[2]], gene %in% 
            rows)
    if (is.null(colors) || "cog.code" %in% names(genome.info)) {
        tmp.lett <- 1:26
        names(tmp.lett) <- LETTERS
        if (!is.null(cluster$cog.code)) 
            coo <- cluster$cog.code[rows]
        else coo <- 1:length(rows)
        tmp <- unique(tmp.lett[coo])
        names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
        cols <- rainbow(length(tmp))
        names(cols) <- names(tmp)
        cols <- cols[names(tmp.lett[coo])]
        cols[is.na(names(cols))] <- "darkgrey"
        names(cols) <- rows
    }
    else {
        cols <- rainbow(length(rows))
        names(cols) <- rows
    }
    cluster$colors <- colors <- cols
    no.motif <- FALSE
    p.values <- motif.widths <- pssm <- NULL
    if (!is.null(motif.out) && !is.null(motif.info) && nrow(motif.info) > 
        0 && (!all(is.na(p.clust)) || !all(is.na(e.clust)))) {
        p.values <- motif.out$p.values[rows]
        motif.widths <- sapply(motif.out$pssms, nrow, simplify = T)
        pssm <- motif.out$pssms
    }
    else {
        no.motif <- TRUE
        p.values <- numeric(length(rows))
        motif.widths <- 0
    }
    seqs <- seqs[rows]
    names(seqs) <- names(p.values) <- rows
    seq.lengths <- nchar(seqs)
    seq.lengths[seq.lengths == 2] <- NA
    if (any(seq.lengths[!is.na(seq.lengths)] > median(seq.lengths, 
        na.rm = T))) {
        seqs <- substr(seqs, 1, median(seq.lengths, na.rm = T))
        seq.lengths <- nchar(seqs)
    }
    maxlen <- max(seq.lengths, na.rm = T)
    if (!is.null(seq.type) && (maxlen == 0 || is.infinite(maxlen))) 
        maxlen <- diff(motif.upstream.search[[seq.type]])
    inds <- integer()
    if (no.motif && (sort.by == "p.value" || sort.by == TRUE)) 
        sort.by <- "gene.name"
    if (sort.by == "gene.name") 
        inds <- sort(rows, decreasing = T, index = T)$ix
    else if (sort.by == "p.value" || sort.by == TRUE) 
        inds <- order(p.values[rows], decreasing = T, na.last = F)
    else if (sort.by == "resid") 
        inds <- order(row.scores[rows, k], decreasing = T)
    if (length(inds) < length(rows)) 
        inds <- c((1:length(rows))[!1:length(rows) %in% inds], 
            inds)
    x.range <- c(-maxlen * 0.08, maxlen * 1.15)
    y.range <- c(0.5, length(rows) + 1)
    plot(x.range, y.range, type = "n", axes = F, xlab = "sequence position", 
        ylab = "")
    cexes <- 1
    if (!no.key) 
        axis(side = 1, pos = 0.6, tck = 0.01, mgp = c(0.1, 0.1, 
            0.1), labels = c(-1, seq(-100, -maxlen, -100)), at = seq(maxlen, 
            0, -100) + motif.upstream.scan[[seq.type]][1], ...)
    if (max(seq.lengths, na.rm = T) > 0) 
        sapply(maxlen - c(0, motif.upstream.search[[seq.type]]) + 
            motif.upstream.scan[[seq.type]][1], function(i) lines(rep(i, 
            2), c(-999, 999), col = "lightgray", lty = 2))
    colmap <- rainbow(length(rows))
    mots.used <- numeric()
    if (is.list(motif.widths)) {
        if (length(motif.widths) <= 0) 
            motif.widths <- 0
        else {
            for (i in 1:length(motif.widths)) if (is.null(motif.widths[[i]])) 
                motif.widths[[i]] <- 0
            motif.widths <- unlist(motif.widths)
        }
    }
    lwd <- 3
    if (length(rows) > 20) 
        lwd <- 1
    else if (length(rows) > 10) 
        lwd <- 2
    if (no.key) 
        lwd <- 1
    if (!no.motif) {
        tmp.mot.info <- subset(motif.info, gene %in% rows)
        tmp.mot.info <- subset(tmp.mot.info, posns <= diff(motif.upstream.scan[[seq.type]]))
        p.min <- quantile(log10(tmp.mot.info$pvals), 0.1, na.rm = T)
        if (is.na(p.min)) 
            p.min <- -5
        p.max <- quantile(log10(tmp.mot.info$pvals), na.rm = T, 
            0.9)
        if (is.na(p.max)) 
            p.max <- log10(p.val.shade.cutoff)
    }
    for (j in 1:length(rows)) {
        jj <- inds[j]
        cur.gene <- rows[jj]
        seq.len <- seq.lengths[jj]
        if (!is.null(rows)) {
            label <- rows[jj]
            if (!is.null(colors)) 
                rect(maxlen + 5, j - 0.18, maxlen * 1.195, j + 
                  0.18, col = colors[label], border = colors[label], 
                  lwd = 3)
        }
        if (!no.motif) {
            rects <- NULL
            mot.info <- subset(tmp.mot.info, gene == cur.gene)
            if (nrow(mot.info) > 0) {
                mots <- mot.info$mots
                starts <- mot.info$posns
                widths <- motif.widths[abs(mots)]
                for (i in 1:length(mots)) {
                  mot <- mots[i]
                  if (is.na(mot) || is.na(seq.len)) 
                    next
                  start <- starts[i]
                  if (start > seq.len) 
                    next
                  end <- start + widths[i]
                  mots.used <- unique(c(mots.used, abs(mot)))
                  col <- abs(mot) + 1
                  if (shade) {
                    if (!is.null(mot.info)) 
                      p.val <- mot.info$pvals[i]
                    else p.val <- 1e-05
                    if (is.na(p.val) || p.val > p.val.shade.cutoff || 
                      p.val > 1) 
                      next
                    else if (p.val <= 0) 
                      p.val <- 1e-05
                    p.val <- log10(p.val)
                    col <- col2rgb(palette()[col])/255
                    col[col > 0] <- 1
                    tmp <- if (p.val < 10) 
                      min(1, max(0, (p.val - p.min)/(p.max - 
                        p.min)))
                    else 0.99
                    if (names(dev.cur()) != "X11") {
                      alpha <- tmp
                    }
                    else {
                      col[col == 0] <- tmp
                      alpha <- 0
                    }
                    col[col < 0] <- 0
                    col[col > 1] <- 1
                    col <- rgb(col["red", 1], col["green", 1], 
                      col["blue", 1], 1 - alpha)
                  }
                  start.1 <- start + maxlen - seq.len
                  end.1 <- end + maxlen - seq.len
                  if (names(dev.cur()) == "devSVG") {
                    par(family = "Arial")
                    setSVGShapeToolTip(title = sprintf("Motif # %2d", 
                      abs(mot)), desc1 = paste(ifelse(mot < 0, 
                      "Rev.", "For."), "strand,", start - maxlen, 
                      "to", end - maxlen), desc2 = sprintf("p-value = %.2g", 
                      10^p.val))
                  }
                  if (!is.null(mot.info)) {
                    if (names(dev.cur()) == "devSVG") {
                      if (mot > 0) 
                        rect(start.1, j + 0.01, end.1, j + 0.3, 
                          col = col, border = col)
                      else if (mot < 0) 
                        rect(start.1, j - 0.3, end.1, j - 0.01, 
                          col = col, border = col)
                    }
                    else {
                      if (mot > 0) 
                        rects <- rbind(rects, c(start.1, j + 
                          0.01, end.1, j + 0.3, col = col, border = col))
                      else if (mot < 0) 
                        rects <- rbind(rects, c(start.1, j - 
                          0.3, end.1, j - 0.01, col = col, border = col))
                    }
                  }
                  else {
                    if (names(dev.cur()) == "devSVG") {
                      if (mot > 0) 
                        rect(start.1, j + 0.01, end.1, j + 0.3, 
                          border = col)
                      else if (mot < 0) 
                        rect(start.1, j - 0.3, end.1, j - 0.01, 
                          border = col)
                    }
                    else {
                      if (mot > 0) 
                        rects <- rbind(rects, c(start.1, j + 
                          0.01, end.1, j + 0.3, col = NA, border = col))
                      else if (mot < 0) 
                        rects <- rbind(rects, c(start.1, j - 
                          0.3, end.1, j - 0.01, col = NA, border = col))
                    }
                  }
                }
            }
            if (!is.null(rects)) 
                rect(rects[, 1], rects[, 2], rects[, 3], rects[, 
                  4], col = rects[, 5], border = rects[, 6])
        }
        slen <- seq.lengths[jj]
        if (all(seq.lengths[!is.na(seq.lengths)]) == 0) 
            slen <- 50
        lines(c(maxlen - slen, maxlen), c(j, j), lwd = lwd + 
            as.integer(rows[jj] %in% o.genes), col = colmap[jj])
        if (grepl("N", seqs[jj])) {
            locs <- which(strsplit(seqs[jj], "")[[1]] == "N")
            diff.locs <- c(diff(locs), 999)
            for (i in 1:sum(diff.locs > 1)) {
                if (i == 1) 
                  lines(c(locs[1], locs[diff.locs > 1][i]) + 
                    maxlen - slen, c(j, j), lwd = lwd + as.integer(rows[jj] %in% 
                    o.genes), col = "gray", lty = 2)
                else lines(c(locs[which(diff.locs > 1) + 1][i - 
                  1], locs[diff.locs > 1][i]) + maxlen - slen, 
                  c(j, j), lwd = lwd + as.integer(rows[jj] %in% 
                    o.genes), col = "gray", lty = 2)
            }
        }
        if (!is.null(rows)) {
            label <- rows[jj]
            col <- "black"
            if (exists("all.tfs") && label %in% all.tfs) 
                col <- "tomato3"
            if (names(dev.cur()) == "devSVG" || !long.names && 
                !no.key) {
                label <- substr(label, 1, 80)
                text(maxlen * 1.2, j, labels = label, adj = c(1, 
                  0.5), col = col, xpd = NA, ...)
            }
            if (long.names || names(dev.cur()) == "devSVG") {
                g.name <- toupper(label)
                if (!is.null(cluster$gene.coords)) 
                  g.name <- cluster$gene.coords[label]
                g.name[is.na(g.name)] <- label[is.na(g.name)]
                if (is.na(g.name)) 
                  g.name <- label
                if (names(dev.cur()) == "devSVG") {
                  par(family = "Arial")
                  setSVGShapeURL(paste("http://www.genome.ad.jp/dbget-bin/www_bget?", 
                    organism, ":", label, sep = ""))
                  if (!is.na(g.name)) 
                    setSVGShapeToolTip(label, g.name)
                  else setSVGShapeToolTip(label)
                  rect(maxlen * 1.2, j - 0.18, maxlen, j + 0.18, 
                    col = NA, border = NA, xpd = NA)
                }
                else if (!is.na(g.name)) {
                  lab <- label
                  if (toupper(g.name) != toupper(label) && g.name != 
                    "") {
                    g.name <- gsub("^[:\\s+]+", "", gsub("\\s+$", 
                      "", g.name, perl = T), perl = T)
                    if (names(dev.cur()) == "X11") 
                      g.name <- strtrim(g.name, 40)
                    else gname <- strtrim(g.name, 60)
                    if (label != "") 
                      lab <- paste(g.name, ": ", label, sep = "")
                    else lab <- g.name
                  }
                  if (nchar(lab) > 60) 
                    lab <- substr(lab, nchar(lab) - 60, nchar(lab))
                  if (!no.key) 
                    text(maxlen * 1.2, j, labels = lab, adj = c(1, 
                      0.5), col = col, xpd = NA, ...)
                }
            }
            if ((!all(is.na(p.clust)) || !all(is.na(e.clust))) && 
                !no.key) 
                text(-maxlen * 0.07, j, labels = sprintf("%.2f", 
                  p.values[label]), xpd = NA, col = if (label %in% 
                  names(which(!is.dup.seq))) 
                  "black"
                else "blue", ...)
        }
    }
    if (!no.key && (!all(is.na(p.clust)) || !all(is.na(e.clust)))) {
        text(-maxlen * 0.15, length(rows) + 0.9, labels = sprintf("log10(P) %s", 
            seq.type), pos = 4, ...)
        mots.used <- sort(unique(mots.used))
        if (length(mots.used) > 1) {
            sapply(1:length(mots.used), function(j) text(maxlen * 
                0.24 + (j + 0) * maxlen * 0.03, length(rows) + 
                0.9, as.character(mots.used[j]), col = mots.used[j] + 
                1, xpd = NA, adj = c(0, 0.5), ...))
        }
        n.unique.seqs <- sum(!is.dup.seq)
        text(maxlen * 1.2, length(rows) + 0.9, sprintf("log10(P.clust)=%.2f; %d seqs; %d uniq", 
            p.clust[seq.type], length(seqs), n.unique.seqs), 
            xpd = NA, adj = c(1, 0.5), ...)
    }
}
plotScores <-
function (k, o.genes = NULL, b.genes = NULL, recompute = F) 
{
    opar <- par(no.readonly = T)
    rows <- get.rows(k)
    if (recompute || !exists("row.scores") || is.null(row.scores)) {
        if (attr(get.all.scores, "version") == 1) {
            tmp <- get.all.scores(k, force.row = T, force.col = T, 
                force.motif = T, force.net = T)
            rs <- tmp$r
            ms <- tmp$m
            ns <- tmp$n
            cs <- tmp$c
        }
        else if (attr(get.all.scores, "version") == 2) {
            tmp <- get.old.scores.matrices(k)
            rs <- tmp$r[, 1]
            ms <- tmp$m[, 1]
            ns <- tmp$n[, 1]
            if (all(is.na(ms)) || all(ms == ms[1])) 
                rm(ms)
            if (all(is.na(ns)) || all(ns == ns[1])) 
                rm(ns)
        }
    }
    tmp.scale <- round(attr(ratios, "nrow")/length(rows)/4)
    layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3, byrow = T))
    if (!exists("rs")) 
        rs <- row.scores[, k, drop = T]
    rs[rs < -220] <- min(rs[rs > -220], na.rm = T)
    h <- try(hist(rs, breaks = 20, main = paste("Cluster", k), 
        xlab = "Ratios scores"), silent = T)
    if (class(h) != "try-error") {
        try(hist(rep(rs[rows], tmp.scale), breaks = h$breaks, 
            col = "red", border = "red", add = T), silent = T)
        try(hist(rs, breaks = h$breaks, add = T), silent = T)
    }
    if (exists("ms") || (!is.null(mot.scores) && !all(is.na(mot.scores[, 
        k])) && !no.genome.info)) {
        if (!exists("ms")) 
            ms <- mot.scores[, k, drop = T]
        ms[ms < -20] <- min(ms[ms > -20], na.rm = T)
        h <- try(hist(ms, breaks = 20, main = NULL, xlab = "Motif scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ms[rows], tmp.scale * 3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ms, breaks = h$breaks, add = T), silent = T)
        }
    }
    else {
        ms <- NULL
        plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
    }
    if (exists("ns") || (!is.null(net.scores) && !all(net.scores[, 
        k] == 0))) {
        if (!exists("ns")) {
            ns <- net.scores[, k, drop = T]
            ns[ns < -20] <- min(ns[ns > -20], na.rm = T)
            ns <- -log10(-ns)
        }
        ns[is.infinite(ns)] <- max(ns[!is.infinite(ns)], na.rm = T) + 
            0.1
        h <- try(hist(ns, breaks = 20, main = NULL, xlab = "-log10(-Network scores)"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ns[rows], tmp.scale/3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ns, breaks = h$breaks, add = T), silent = T)
        }
    }
    else {
        ns <= NULL
        plot(1, 1, typ = "n", axes = F, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
    }
    row.memb <- attr(ratios, "rnames") %in% rows
    if (!is.null(ms) && !all(is.na(ms)) && !all(ns == 0)) {
        plot(rs, ms, typ = "n", main = paste("Cluster", k), xlab = "Ratios scores", 
            ylab = "Mot scores")
        text(rs, ms, label = 1:length(rs), col = row.memb + 1, 
            cex = 0.5)
    }
    else if (!is.null(ns) && !all(ns == 0)) {
        plot(rs, ns, typ = "n", main = paste("Cluster", k), xlab = "Ratios scores", 
            ylab = "Net scores")
        text(rs, ns, label = 1:length(rs), col = row.memb + 1, 
            cex = 0.5)
    }
    else {
        plot(rs, jitter(rep(0, length(rs))), typ = "n", main = paste("Cluster", 
            k), xlab = "Ratios scores", ylab = "")
        text(rs, jitter(rep(0, length(rs))), label = 1:length(rs), 
            col = row.memb + 1, cex = 0.5)
    }
    if (!is.null(o.genes)) 
        text(rs[o.genes], ms[o.genes], label = which(attr(ratios, 
            "rnames") %in% o.genes), col = "green", cex = 0.5)
    if (!is.null(b.genes)) 
        text(rs[b.genes], ms[b.genes], label = which(attr(ratios, 
            "rnames") %in% b.genes), col = "blue", cex = 0.5)
    try({
        tmp <- get.combined.scores(quant = T)
        r.scores <- tmp$r
        c.scores <- tmp$c
        rr <- get.density.scores(ks = k, r.scores, c.scores, 
            plot = "rows")$r
        rr <- rr[, k, drop = T]
        h <- try(hist(log10(rr), breaks = 50, main = NULL, xlab = "Density (membership) scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(log10(rr[rows]), tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(log10(rr), breaks = h$breaks, add = T), 
                silent = T)
        }
    }, silent = T)
    par(opar)
}
plotStats <-
function (iter = stats$iter[nrow(stats)], plot.clust = NA, new.dev = F, 
    ...) 
{
    if (!exists("row.memb")) {
        row.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "rnames") %in% get.rows(k))
        if (is.vector(row.memb)) 
            row.memb <- t(row.memb)
        rownames(row.memb) <- attr(ratios, "rnames")
        col.memb <- sapply(1:k.clust, function(k) attr(ratios, 
            "cnames") %in% get.cols(k))
        if (is.vector(col.memb)) 
            col.memb <- t(col.memb)
        rownames(col.memb) <- attr(ratios, "cnames")
    }
    if (!exists("row.scores") || is.null(row.scores)) {
        if (attr(get.all.scores, "version") == 1) {
            tmp <- get.all.scores()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
            tmp <- get.combined.scores(quant = T)
            r.scores <- tmp$r
            c.scores <- tmp$c
        }
        else if (attr(get.all.scores, "version") == 2) {
            tmp <- get.old.scores.matrices()
            row.scores <- tmp$r
            mot.scores <- tmp$m
            net.scores <- tmp$n
            col.scores <- tmp$c
        }
    }
    opar <- par(no.readonly = T)
    tmp.scale <- round(1/mean(row.memb, na.rm = T)/4)
    if (new.dev) {
        if (length(dev.list()) < 1) 
            dev.new()
        dev.set(2)
    }
    layout(matrix(c(1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 7, 9, 
        11, 8, 10, 12), byrow = T, ncol = 3))
    par(mar = c(3, 3, 2, 0.1), mgp = c(3, 1, 0) * 0.5)
    stats <- stats[stats[, "iter"] <= iter, , drop = F]
    try(matplot(stats[, "iter"], stats[, grep("resid", colnames(stats), 
        val = T)], typ = "l", xlab = "iter", ylab = "Mean resid", 
        main = sprintf("Iter: %d", iter), lty = 1), silent = T)
    sapply(c(51, 101, 21), function(kwin) try(matlines(stats[, 
        "iter"], apply(stats[, grep("resid", colnames(stats), 
        val = T), drop = F], 2, function(i) runmed(i, k = min(length(i), 
        kwin))), lty = 2, lwd = 0.6), silent = T))
    if ((nn <- length(grep("resid", colnames(stats)))) > 1) 
        legend("bottomleft", legend = gsub("resid.", "", grep("resid", 
            colnames(stats), val = T)), lwd = 1, bty = "n", col = 1:nn, 
            lty = 1:nn, cex = 0.5)
    if (exists("row.scores") && !is.null(mot.scores) && !all(is.na(mot.scores[, 
        ]))) {
        rs <- row.scores[]
        rs[rs < -20] <- min(rs[rs > -20], na.rm = T)
        h <- try(hist(rs, breaks = 50, main = NULL, xlab = "Ratios scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(rs[row.memb == 1], tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(rs, breaks = h$breaks, add = T), silent = T)
        }
    }
    if (exists("mot.scores") && !is.null(mot.scores) && !all(is.na(mot.scores[, 
        ]))) {
        ms <- mot.scores[, ]
        ms[ms < -20] <- min(ms[ms > -20], na.rm = T)
        ms[ms >= 0] <- NA
        h <- try(hist(ms, breaks = 50, main = NULL, xlab = "Motif scores"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ms[row.memb == 1], tmp.scale * 3), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ms, breaks = h$breaks, add = T), silent = T)
        }
    }
    tmp <- stats[, grep("p.clust", colnames(stats), val = T), 
        drop = F]
    if (!all(is.na(tmp))) {
        try(matplot(stats[, "iter"], tmp, typ = "l", xlab = "iter", 
            ylab = "Mean motif p-value", main = sprintf("Motif scaling: %.3f", 
                mot.scaling[max(1, iter - 1)]), lty = 1), silent = T)
        sapply(c(51, 101, 21), function(kwin) try(matlines(stats[!is.na(tmp), 
            "iter"], apply(tmp, 2, function(i) runmed(i[!is.na(i)], 
            k = min(sum(!is.na(i)), kwin))), lty = 2, lwd = 0.6), 
            silent = T))
        if ((nn <- length(grep("p.clust", colnames(stats)))) > 
            1) 
            legend("bottomleft", legend = gsub("p.clust.", "", 
                grep("p.clust", colnames(stats), val = T)), lwd = 1, 
                bty = "n", col = 1:nn, lty = 1:nn, cex = 0.5)
    }
    if (exists("net.scores") && !is.null(net.scores) && !all(net.scores[, 
        ] == 0)) {
        ns <- net.scores[, ]
        ns[ns < -20] <- min(ns[ns > -20], na.rm = T)
        ns[ns >= 0] <- NA
        ns[, ] <- -log10(-ns)
        tmp.scale <- ceiling(tmp.scale * mean(!is.na(ns), na.rm = T))
        h <- try(hist(ns, breaks = 50, main = NULL, xlab = "-log10(-Network scores)"), 
            silent = T)
        if (class(h) != "try-error") {
            try(hist(rep(ns[row.memb == 1], tmp.scale), breaks = h$breaks, 
                col = "red", border = "red", add = T), silent = T)
            try(hist(ns, breaks = h$breaks, add = T), silent = T)
        }
        tmp <- stats[, grep("net.", colnames(stats), val = T, 
            fixed = T), drop = F]
        try(matplot(stats[, "iter"], tmp, typ = "l", xlab = "iter", 
            ylab = "Mean net-score", main = sprintf("Net scaling: %.3f", 
                net.scaling[max(1, iter - 1)]), lty = 1), silent = T)
        try(matlines(stats[, "iter"], apply(tmp, 2, function(i) runmed(i[!is.na(i)], 
            k = min(sum(!is.na(i)), 51))), lty = 2, lwd = 0.6), 
            silent = T)
        if ((nn <- length(grep("net.", colnames(stats)))) > 1) 
            try(legend("bottomleft", legend = gsub("net.", "", 
                grep("net.", colnames(stats), val = T)), lwd = 1, 
                bty = "n", col = 1:nn, lty = 1:nn, cex = 0.5), 
                silent = T)
    }
    clusterStack <- get.clusterStack(ks = 1:k.clust)
    resids <- sapply(as.list(clusterStack), "[[", "resid")
    try(hist(resids[resids <= 1.5], main = NULL, xlab = "Cluster Residuals", 
        xlim = if (all(resids > 0, na.rm = T)) 
            c(0, 1.5)
        else range(resids, na.rm = T), breaks = k.clust/4), silent = T)
    if (!exists("mot.scores") || is.null(mot.scores) || all(is.na(mot.scores[, 
        ]))) {
        pclusts <- sapply(as.list(clusterStack), "[[", "p.clust")
        try(hist(pclusts[pclusts <= 1], main = NULL, xlab = "Cluster Motif P-values", 
            xlim = if (all(pclusts > 0, na.rm = T)) 
                c(0, 1)
            else range(pclusts, na.rm = T), breaks = k.clust/4), 
            silent = T)
    }
    if (mot.scaling[iter] > 0) {
        plot.all.clusterMotifPositions <- function(ks = 1:k.clust, 
            mots = 1, e.cutoff = 1, p.cutoff = 0.05, seq.type = names(mot.weights)[1], 
            breaks = 100, ...) {
            if (seq.type == "ALL") 
                seq.type <- names(mot.weights)
            df <- NULL
            for (st in seq.type) {
                ms <- meme.scores[[st]][ks]
                ind <- 1
                if (!"posns" %in% colnames(ms[[1]]$pv.ev[[1]])) 
                  ind <- 2
                posns <- as.vector(unlist(sapply(ms, function(i) i$pv.ev[[ind]]$posns)))
                pvals <- as.vector(unlist(sapply(ms, function(i) i$pv.ev[[ind]]$pvals)))
                imots <- as.vector(unlist(sapply(ms, function(i) i$pv.ev[[ind]]$mots)))
                genes <- as.vector(unlist(sapply(ms, function(i) as.character(i$pv.ev[[ind]]$gene))))
                slens <- nchar(genome.info$all.upstream.seqs[st][[st]][genes])
                clusts <- as.vector(unlist(sapply(ms, function(i) rep(i$k, 
                  if (is.null(i$pv.ev[[ind]])) 0 else nrow(i$pv.ev[[ind]])))))
                ms <- meme.scores[[st]]
                evals <- sapply(1:length(imots), function(i) ms[[clusts[i]]]$meme.out[[abs(imots[i])]]$e.value)
                df <- rbind(df, data.frame(clusts, posns, pvals, 
                  imots, evals, genes, slens, seq.type = rep(st, 
                    length(clusts))))
            }
            df2 <- subset(df, evals < e.cutoff & pvals < p.cutoff & 
                abs(imots) %in% mots)
            psns <- df2$posns - df2$slens - do.call(rbind, motif.upstream.scan[df2$seq.type])[, 
                1]
            h <- hist(psns, breaks = breaks, xlab = sprintf("%s %s", 
                paste(seq.type, collapse = " "), paste(mots, 
                  collapse = " ")), ...)
            dd <- density(psns, bw = 5)
            lines(dd$x, dd$y * max(h$counts)/max(dd$y) * 0.9, 
                col = "red")
            invisible(df)
        }
        try(plot.all.clusterMotifPositions(main = "Positions of motif #1", 
            ...), silent = T)
        if (any(sapply(meme.scores$upstream[1:k.clust], function(i) length(i$meme.out)) == 
            2)) 
            try(plot.all.clusterMotifPositions(mots = 2, main = "Positions of motif #2", 
                ...), silent = T)
    }
    n.rows <- sapply(1:k.clust, function(k) length(get.rows(k)))
    try(hist(n.rows, main = NULL, xlab = "Cluster Nrows", breaks = k.clust/4, 
        xlim = c(-5, max(n.rows, na.rm = T))), silent = T)
    n.cols <- sapply(1:k.clust, function(k) length(get.cols(k)))
    try(hist(n.cols, main = NULL, xlab = "Cluster Ncols", breaks = k.clust/4, 
        xlim = c(-5, attr(ratios, "ncol"))), silent = T)
    nr <- table(unlist(sapply(clusterStack, "[[", "rows")))
    if (length(nr) < attr(ratios, "nrow")) 
        nr <- c(nr, rep(0, attr(ratios, "nrow") - length(nr) + 
            1))
    try(hist(nr, breaks = seq(-0.5, 10, by = 1), main = NULL, 
        xlab = "NClust per gene"), silent = T)
    if (!is.na(plot.clust)) {
        if (new.dev) {
            if (length(dev.list()) < 2) 
                dev.new()
            dev.set(3)
        }
        try(plotClust(plot.clust, w.motifs = T, cex = 0.7), silent = T)
        if (new.dev) {
            if (length(dev.list()) < 3) 
                dev.new()
            dev.set(4)
        }
        try(plotScores(plot.clust), silent = T)
    }
    par(opar)
}
preprocess.ratios <-
function (ratios, filter = T, normalize = T, col.groups = NULL, 
    frac.cutoff = 0.98) 
{
    if (is.null(col.groups)) 
        col.groups <- rep(1, ncol(ratios))
    if (is.null(names(col.groups))) 
        names(col.groups) <- colnames(ratios)
    if (filter) {
        cat("Filtering out nochange rows/cols from ratios matrix...\n")
        tmp1 <- apply(ratios, 1, function(i) mean(is.na(i) | 
            abs(i) <= 0.17)) < frac.cutoff
        tmp2 <- apply(ratios, 2, function(i) mean(is.na(i) | 
            abs(i) <= 0.1)) < frac.cutoff
        ratios <- ratios[tmp1, , drop = F]
        ratios <- ratios[, tmp2, drop = F]
        cat("Filtered ratios matrix is", paste(dim(ratios), collapse = "x"), 
            "\n")
        col.groups <- col.groups[tmp2]
    }
    if (normalize) {
        for (cg in unique(col.groups)) {
            cols <- names(which(col.groups == cg))
            cat("Converting ratios matrix", cg, "to z-scores...\n")
            ratios[, cols] <- t(scale(t(ratios[, cols, drop = F]), 
                center = apply(ratios[, cols, drop = F], 1, median, 
                  na.rm = T), scale = apply(ratios[, cols, drop = F], 
                  1, sd, na.rm = T)))
        }
    }
    ratios
}
pssm.motif.lines <-
function (pssm, id, e.value = 1, header = T, seq.type = "upstream weeder") 
{
    meme.let <- c("A", "C", "G", "T")
    if (missing(id)) 
        id <- paste(pssm, collapse = "")
    lines <- character()
    if (header) {
        lines <- "ALPHABET= ACGT"
        lines <- c(lines, paste(names(unlist(genome.info$bg.list[[seq.type]][meme.let])), 
            sprintf("%.3f", unlist(genome.info$bg.list[[seq.type]][meme.let])), 
            collapse = " "))
    }
    if (is.null(colnames(pssm))) 
        colnames(pssm) <- col.let
    pssm <- pssm[, meme.let] + max(pssm, na.rm = T)/100
    for (i in 1:nrow(pssm)) pssm[i, ] <- pssm[i, ]/sum(pssm[i, 
        ], na.rm = T)
    idd <- gsub("[_/]", ".", id)
    lines <- c(lines, sprintf("log-odds matrix: alength= 4 w= %d", 
        nrow(pssm)))
    for (j in 1:nrow(pssm)) lines <- c(lines, paste(sprintf("%5.3f", 
        log2(pssm[j, ])), collapse = " ", sep = " "))
    lines
}
pssm.to.consensus <-
function (pssm, cutoff.1 = 0.8, cutoff.2 = 0.6, regex = T) 
{
    c1 <- apply(pssm, 1, function(i) which(i > cutoff.1))
    c2 <- apply(pssm, 1, function(i) which(i > cutoff.2))
    ca <- rbind(sapply(c1, length), sapply(c2, length))
    cond <- function(l) paste("[", paste(l, collapse = ""), "]", 
        sep = "")
    deg.codes <- c(AG = "R", GT = "K", CG = "S", CT = "Y", AC = "M", 
        AT = "W", CGT = "B", ACT = "H", AGT = "D", ACG = "V", 
        ACGT = "N")
    out <- character()
    for (i in 1:length(c1)) {
        if (ca[1, i] == 1) {
            out[i] <- col.let[c1[[i]]]
        }
        else if (ca[2, i] >= 2) {
            if (sum(pssm[i, c2[[i]]]) > cutoff.1) 
                out[i] <- cond(col.let[c2[[i]]])
            else if (sum(pssm[i, c2[[i]]]) > cutoff.2) 
                out[i] <- cond(tolower(col.let[c2[[i]]]))
        }
        else if (ca[2, i] >= 1) {
            if (pssm[i, c2[[i]]] > cutoff.2) 
                out[i] <- tolower(col.let[[c2[[i]]]])
        }
        if (is.na(out[i])) 
            out[i] <- "[ACGT]"
        if (!regex && nchar(out[i]) > 1) {
            tmp <- substr(out[i], 2, nchar(out[i]) - 1)
            if (tmp %in% names(deg.codes)) 
                out[i] <- deg.codes[tmp]
            else if (toupper(tmp) %in% names(deg.codes)) 
                out[i] <- tolower(deg.codes[toupper(tmp)])
        }
    }
    paste(out, collapse = "")
}
pssm.to.string <-
function (pssm, cutoff.1 = 0.7, cutoff.2 = 0.4) 
{
    maxes <- max.col(pssm)
    letters <- col.let[maxes]
    values <- pssm[cbind(1:nrow(pssm), maxes)]
    letters[letters == "A" & values < cutoff.1] <- "a"
    letters[letters == "C" & values < cutoff.1] <- "c"
    letters[letters == "G" & values < cutoff.1] <- "g"
    letters[letters == "T" & values < cutoff.1] <- "t"
    letters[values < cutoff.2] <- "n"
    return(paste(letters, collapse = ""))
}
quantile.normalize.scores <-
function (scores, weights = NULL, keep.nas = F) 
{
    if (!is.list(scores) || sum(!sapply(scores, is.null)) <= 
        1) 
        return(scores)
    scores <- scores[!sapply(scores, is.null)]
    d <- dim(scores[[1]])
    dn <- dimnames(scores[[1]])
    tmp2 <- sapply(scores, function(i) sort(i[, ], na.last = T))
    if (is.null(weights)) 
        tmp2.mn <- rowMeans(tmp2, na.rm = T)
    else {
        for (i in 1:ncol(tmp2)) tmp2[, i] <- tmp2[, i] * weights[i]
        tmp2.mn <- rowMeans(tmp2, na.rm = T)/sum(weights, na.rm = T)
    }
    rm(tmp2)
    out <- list()
    for (n in names(scores)) {
        tmp <- rank(scores[[n]][, ], ties = "min", na = "keep")
        z <- matrix(tmp2.mn[tmp], nrow = d[1], ncol = d[2])
        dimnames(z) <- dn
        if (keep.nas) 
            z[is.na(scores[[n]][, ])] <- NA
        out[[n]] <- z
    }
    out
}
re.seed.empty.clusters <-
function (row.membership, col.membership, toosmall.r = cluster.rows.allowed[1], 
    toosmall.c = 0, toobig.r = cluster.rows.allowed[2], n.r = cluster.rows.allowed[1] * 
        2, n.c = 5) 
{
    rm <- row.membership
    rats <- get.cluster.matrix()
    if (any(tabulate(unlist(apply(rm, 1, unique)), k.clust) <= 
        toosmall.r)) {
        which.zero <- which(tabulate(unlist(apply(rm, 1, unique)), 
            k.clust) <= toosmall.r)
        cat("These", length(which.zero), "clusters have TOO FEW rows: ", 
            which.zero, "\n")
        which.toobig <- which(tabulate(unlist(apply(rm, 1, unique)), 
            k.clust) >= toobig.r)
        cat("These", length(which.toobig), "clusters have TOO MANY rows: ", 
            which.toobig, "\n")
        which.zero <- c(which.zero, which.toobig)
        for (k in which.zero) {
            all.zero <- names(which(apply(rm, 1, function(i) all(i <= 
                toosmall.r))))
            if (length(all.zero) < n.r) {
                all.zero <- unique(c(all.zero, rownames(which(rm == 
                  0, arr = T))))
                all.zero <- unique(c(all.zero, names(which(apply(rm, 
                  1, function(i) all(i == i[1]))))))
            }
            if (length(all.zero) <= 1) 
                break
            gs <- sample(all.zero, 1)
            cors <- apply(rats[all.zero, ], 1, cor, rats[gs, 
                ], use = "pairwise")
            gs <- names(cors[order(cors, decreasing = T)[1:n.r]])
            gs <- gs[!is.na(gs)]
            for (g in gs) {
                if (any(rm[g, ] == 0)) 
                  rm[g, which(rm[g, ] == 0)[1]] <- k
                else rm[g, 1] <- k
            }
        }
        for (tt in names(mot.weights)) for (k in which.zero) meme.scores[[tt]][[k]] <- list(iter = iter)
    }
    cm <- col.membership
    if (any(tabulate(cm, k.clust) <= toosmall.c)) {
        which.zero <- which(tabulate(cm, k.clust) <= toosmall.c)
        cat("These", length(which.zero), "clusters have TOO FEW columns: ", 
            which.zero, "\n")
        for (k in which.zero) {
            all.zero <- names(which(apply(cm, 1, function(i) all(i <= 
                toosmall.c))))
            if (length(all.zero) <= n.c) 
                all.zero <- unique(c(all.zero, rownames(which(cm == 
                  0, arr = T))))
            if (length(all.zero) <= 1) 
                break
            cs <- unique(sample(all.zero, min(length(all.zero), 
                n.c)))
            cs <- cs[!is.na(cs)]
            for (cc in cs) cm[cc, which(cm[cc, ] == 0)[1]] <- k
        }
    }
    invisible(list(r = rm, c = cm, ms = meme.scores))
}
read.fasta <-
function (fname, lines = NULL) 
{
    if (is.null(lines)) 
        lines <- readLines(fname)
    lines <- lines[lines != ""]
    starts <- grep("^>", lines, perl = T)
    if (length(starts) > 1) 
        stops <- c(starts[2:length(starts)], length(lines) + 
            1)
    else stops <- length(lines) + 1
    seqs <- sapply(1:length(starts), function(i) paste(lines[(starts[i] + 
        1):(stops[i] - 1)], collapse = "", sep = ""))
    names(seqs) <- gsub("^>", "", lines[starts], perl = T)
    seqs
}
remove.low.complexity <-
function (seqs, length = 8, entropy.cutoff = 0.6, repl = "N", 
    use.dust = T, seq.type = names(mot.weights)[1]) 
{
    write.fasta <- function(seqs, fname) writeLines(paste(paste(">", 
        names(seqs), sep = ""), seqs, sep = "\n"), con = fname)
    if (use.dust) {
        seqs <- seqs[!is.null(seqs) & !is.na(seqs)]
        max.width <- as.integer(strsplit(meme.cmd[seq.type], 
            " ")[[1]][which(strsplit(meme.cmd[seq.type], " ")[[1]] == 
            "-maxw") + 1])
        seqs <- seqs[nchar(seqs) >= max.width]
        if (length(seqs) > 0) {
            fname <- my.tempfile("dust", suf = ".fst")
            write.fasta(seqs, fname)
            cmd <- gsub("$fname", fname, dust.cmd, fixed = T)
            fst <- system.time.limit(paste(cmd, "2>/dev/null"), 
                tlimit = 60)
            unlink(fname)
            if (length(fst) <= 1) 
                cat("WARNING: you probably don't have 'dust' installed.\n")
            else seqs <- read.fasta(NULL, fst)
            return(seqs)
        }
    }
    shannon.entropy <- function(string) {
        ni <- table(string)/length + 1e-10
        -sum(ni * log2(ni))
    }
    all.dna.seqs <- function(l, lett = c("G", "A", "T", "C"), 
        as.matrix = F) {
        n.lett <- length(lett)
        out <- sapply(1:l, function(ll) rep(as.vector(sapply(lett, 
            function(i) rep(i, n.lett^(ll - 1)))), n.lett^(l - 
            ll)))
        if (as.matrix) 
            return(out)
        apply(out, 1, paste, collapse = "")
    }
    substrings <- all.dna.seqs(l = length)
    mc <- get.parallel(length(substrings))
    bad.substrings <- substrings[unlist(mc$apply(strsplit(substrings, 
        NULL), shannon.entropy)) <= entropy.cutoff]
    repl <- paste(rep(repl, length), sep = "", collapse = "")
    in.seqs <- seqs
    mc <- get.parallel(length(seqs))
    seqs <- unlist(mc$apply(1:length(seqs), function(i) {
        seq <- seqs[i]
        for (s in bad.substrings) seq <- gsub(s, repl, seq)
        seq
    }))
    names(seqs) <- names(in.seqs)
    return(seqs)
}
rev.comp <-
function (seqs) 
{
    sapply(seqs, function(seq) paste(rev(strsplit(toupper(chartr("ATCG", 
        "tagc", seq)), "")[[1]]), collapse = ""))
}
row.col.membership.from.clusterStack <-
function (cs) 
{
    row.memb <- col.memb <- NULL
    for (k in 1:length(cs)) {
        row.memb <- cbind(row.memb, rep(0, attr(ratios, "nrow")))
        if (ncol(row.memb) == 1) 
            rownames(row.memb) <- attr(ratios, "rnames")
        rows <- cs[[k]]$rows
        rows <- rows[!is.na(rows)]
        row.memb[rows, k] <- k
        col.memb <- cbind(col.memb, rep(0, attr(ratios, "ncol")))
        if (ncol(col.memb) == 1) 
            rownames(col.memb) <- attr(ratios, "cnames")
        cols <- cs[[k]]$cols
        cols <- cols[!is.na(cols)]
        col.memb[cols, k] <- k
    }
    row.memb <- t(apply(row.memb, 1, function(i) c(i[i != 0], 
        i[i == 0])))
    row.memb <- row.memb[, apply(row.memb, 2, sum) != 0, drop = F]
    colnames(row.memb) <- NULL
    col.memb <- t(apply(col.memb, 1, function(i) c(i[i != 0], 
        i[i == 0])))
    col.memb <- col.memb[, apply(col.memb, 2, sum) != 0, drop = F]
    colnames(col.memb) <- NULL
    if (ncol(row.memb) < n.clust.per.row) 
        row.memb <- cbind(row.memb, rep(0, nrow(row.memb)))
    if (ncol(col.memb) < n.clust.per.col) 
        col.memb <- cbind(col.memb, rep(0, nrow(col.memb)))
    list(r = row.memb, c = col.memb)
}
runMast <-
function (memeOut, mast.cmd, genes, seqs, bgseqs = NULL, bg.list = NULL, 
    bgfname = NULL, unlink = T, verbose = F, ...) 
{
    fname <- my.tempfile("mast.tmp", suf = ".fst")
    if (is.null(bgfname)) {
        bgfname <- my.tempfile("mast.tmp", suf = ".bg")
        on.exit(unlink(bgfname))
    }
    memeOutFname <- my.tempfile("meme.tmp", suf = ".out")
    cat(memeOut, sep = "\n", file = memeOutFname)
    tmp <- mkTempMemeFiles(genes, seqs, fname = fname, bgseqs = bgseqs, 
        bg.list = bg.list, bgfname = bgfname, seq.weights = NULL, 
        psps = NULL, ...)
    if (tmp <= 0) 
        return(NULL)
    cmd <- mast.cmd
    if (is.null(bgfname) || !file.exists(bgfname)) 
        cmd <- gsub("-bfile $bgFname", "", cmd, fixed = T)
    else cmd <- gsub("$bgFname", bgfname, cmd, fixed = T)
    cmd <- gsub("$memeOutFname", memeOutFname, cmd, fixed = T)
    cmd <- gsub("$fname", fname, cmd, fixed = T)
    if (verbose) 
        cat(cmd, "\n")
    output <- system.time.limit(cmd)
    attr(output, "mast.command.line") <- cmd
    if (unlink) 
        unlink(c(memeOutFname, fname))
    output
}
runMeme <-
function (sgenes, seqs, cmd = meme.cmd[names(mot.weights)[1]], 
    bgseqs = NULL, bgfname = NULL, bg.list = NULL, nmotif = 1, 
    unlink = T, verbose = T, seq.weights = NULL, psps = NULL, 
    ...) 
{
    fname <- my.tempfile("meme.tmp", suf = ".fst")
    if (is.null(bgfname)) {
        bgfname <- my.tempfile("meme.tmp", suf = ".bg")
        on.exit(unlink(bgfname))
    }
    tmp <- mkTempMemeFiles(sgenes, seqs, fname = fname, bgseqs = bgseqs, 
        bg.list = bg.list, bgfname = bgfname, seq.weights = seq.weights, 
        psps = psps, ...)
    if (tmp <= 0) 
        return(NULL)
    if (is.null(bgfname) || !file.exists(bgfname)) 
        cmd <- gsub("-bfile $bgFname", "", cmd, fixed = T)
    else cmd <- gsub("$bgFname", bgfname, cmd, fixed = T)
    if (is.null(psps)) 
        cmd <- gsub("-psp $pspFname", "", cmd, fixed = T)
    else cmd <- gsub("$pspFname", sprintf("%s.psp", fname), cmd, 
        fixed = T)
    cmd <- gsub("$fname", fname, cmd, fixed = T)
    if (verbose) 
        cat(cmd, "\n")
    output <- system.time.limit(cmd)
    attr(output, "meme.command.line") <- cmd
    if (unlink) 
        unlink(c(fname, sprintf("%s.psp", fname)))
    return(output)
}
save.cmonkey.env <-
function (env = NULL, file = NULL, verbose = T) 
{
    if (is.null(env)) {
        for (i in ls(env = globalenv())) if (is.environment(get(i, 
            globalenv())) && "cmonkey" %in% class(get(i, globalenv()))) 
            save.cmonkey.env(get(i, globalenv()), file, verbose)
        return(invisible())
    }
    if (is.null(file)) 
        file <- paste(env$cmonkey.filename, ".RData", sep = "")
    un.ffify.env(env)
    if (verbose) 
        message("Saving environment to ", file)
    save(env, file = file)
    invisible(env)
}
seed.clusters <-
function (k.clust, seed.method = "rnd", col.method = "rnd") 
{
    if (seed.method == "custom" && exists("seed.clusters.custom")) 
        return(seed.clusters.custom(k.clust, col.method))
    if (substr(seed.method, 1, 3) == "net" && length(networks) <= 
        0) {
        cat("Seed method is", seed.method, ", but no networks -- using 'kmeans' instead.\n")
        seed.method <- "kmeans"
    }
    if (seed.method == "rnd") {
        rm <- t(sapply(1:attr(ratios, "nrow"), function(i) sample(1:k.clust, 
            n.clust.per.row[1], replace = n.clust.per.row[1] > 
                attr(ratios, "nrow"))))
    }
    else if (substr(seed.method, 1, 5) == "list=") {
        rm <- matrix(0, nrow = attr(ratios, "nrow"), ncol = n.clust.per.row[1])
        rownames(rm) <- attr(ratios, "rnames")
        fname <- strsplit(seed.method, "=")[[1]][2]
        if (exists(fname)) 
            lists <- get(fname)
        else if (file.exists(fname)) 
            lists <- strsplit(readLines(fname), split = "[ \t,]", 
                perl = T)
        for (k in 1:min(c(k.clust, length(lists)))) {
            probes <- unlist(lapply(get.synonyms(lists[[k]]), 
                function(i) i %in% rownames(rm)))
            rm[probes[rm[probes, 1] == 0], 1] <- k
            rm[probes[rm[probes, 1] != 0], 2] <- k
        }
        if (length(lists) < k.clust) {
            for (k in (length(lists) + 1):k.clust) {
                rnames <- attr(ratios, "rnames")[!attr(ratios, 
                  "rnames") %in% unlist(lists)]
                rows <- sample(rnames, 5)
                rm[rows[rm[rows, 1] == 0], 1] <- k
                rm[rows[rm[rows, 1] != 0 & rm[rows, 2] == 0], 
                  2] <- k
            }
        }
    }
    else if (substr(seed.method, 1, 4) == "rnd=") {
        n.samp <- as.integer(strsplit(seed.method, "=")[[1]][2])
        rm <- matrix(0, nrow = attr(ratios, "nrow"), ncol = n.clust.per.row[1])
        rownames(rm) <- attr(ratios, "rnames")
        for (i in 1:n.clust.per.row) {
            sampled <- rep(FALSE, attr(ratios, "nrow"))
            names(sampled) <- attr(ratios, "rnames")
            for (k in 1:k.clust) {
                g <- sample(attr(ratios, "rnames")[!sampled], 
                  n.samp)
                rm[g, 1] <- k
                sampled[g] <- TRUE
            }
        }
    }
    else if (seed.method == "kmeans") {
        if (!exists("ratios")) 
            stop("kmeans seed method but no ratios")
        tmp.rat <- get.cluster.matrix()
        tmp.rat[is.na(tmp.rat)] <- 0
        km <- kmeans
        if (!is.na(parallel.cores) && (parallel.cores == TRUE || 
            parallel.cores > 1)) {
            get.parallel()
            if (require(biganalytics)) 
                km <- bigkmeans
        }
        rm <- km(tmp.rat, centers = k.clust, iter.max = 20, nstart = 2)$cluster
        names(rm) <- attr(ratios, "rnames")
        if (n.clust.per.row[1] > 1) 
            rm <- cbind(rm, matrix(rep(0, attr(ratios, "nrow") * 
                (n.clust.per.row[1] - 1)), ncol = n.clust.per.row[1] - 
                1))
    }
    else if (substr(seed.method, 1, 11) == "trimkmeans=") {
        if (!exists("ratios")) 
            stop("trimkmeans seed method but no ratios")
        require(trimcluster)
        trim <- as.numeric(strsplit(seed.method, "=")[[1]][2])
        tmp.rat <- get.cluster.matrix()
        tmp.rat[is.na(tmp.rat)] <- 0
        rm <- trimkmeans(tmp.rat, k.clust, trim = trim, maxit = 20, 
            runs = 2)$classification
        if (n.clust.per.row[1] > 1) 
            rm <- cbind(rm, matrix(rep(0, attr(ratios, "nrow") * 
                (n.clust.per.row[1] - 1), ncol = n.clust.per.row[1] - 
                1)))
    }
    else if (substr(seed.method, 1, 4) == "cor=") {
        if (!exists("ratios")) 
            stop("cor seed method but no ratios")
        n.cor <- as.integer(strsplit(seed.method, "=")[[1]][2])
        rats <- get.cluster.matrix()
        cors <- if (attr(ratios, "nrow") < 6000) 
            cor(t(rats), use = "pairwise")
        else NULL
        rm <- rep(0, attr(ratios, "nrow"))
        names(rm) <- attr(ratios, "rnames")
        sampled <- rep(FALSE, attr(ratios, "nrow"))
        names(sampled) <- attr(ratios, "rnames")
        mc <- get.parallel(n.clust.per.row)
        tmp <- mc$apply(1:n.clust.per.row, function(i) {
            for (k in 1:k.clust) {
                if (sum(!sampled) < n.cor) 
                  sampled[sample(1:length(sampled))] <- FALSE
                rnames <- attr(ratios, "rnames")[!sampled]
                g <- sample(rnames, 1)
                if (!is.null(cors)) 
                  g <- rnames[order(cors[g, !sampled], decreasing = T)[1:n.cor]]
                else g <- rnames[order(apply(rats[!sampled, ], 
                  1, cor, rats[g, ]), decreasing = T)[1:n.cor]]
                rm[g] <- k
                if (length(g) == 1) {
                  if (!is.null(cors)) 
                    tmp <- rnames[order(cors[g, !sampled], decreasing = T)[1:10]]
                  else tmp <- rnames[order(apply(rats[!sampled, 
                    ], 1, cor, rats[g, ]), decreasing = T)[1:10]]
                  sampled[tmp] <- TRUE
                }
                sampled[g] <- TRUE
            }
            rm[attr(ratios, "rnames")]
        })
        rm <- do.call(cbind, tmp)
    }
    else if (substr(seed.method, 1, 4) == "net=") {
        if (!exists("networks") || length(networks) <= 0) 
            stop("net seed method but no networks")
        seed.method <- strsplit(seed.method, "=")[[1]][2]
        net.name <- strsplit(seed.method, ":")[[1]][1]
        net <- networks[[net.name]]
        n.seed <- as.integer(strsplit(seed.method, ":")[[1]][2])
        rm <- rep(0, attr(ratios, "nrow"))
        names(rm) <- attr(ratios, "rnames")
        sampled <- rep(FALSE, length(unique(as.character(net$protein1))))
        names(sampled) <- unique(as.character(net$protein1))
        mc <- get.parallel(n.clust.per.row)
        tmp <- mc$apply(1:n.clust.per.row, function(i) {
            for (k in 1:k.clust) {
                if (sum(!sampled) <= 0) 
                  sampled[1:length(sampled)] <- FALSE
                rnames <- names(which(!sampled))
                gs <- sample(rnames, 1)
                qiter <- 0
                while (length(gs) < n.seed && qiter < 20) {
                  ns <- as.character(net$protein2[as.character(net$protein1) %in% 
                    gs])
                  if (length(ns) + length(gs) >= n.seed) 
                    ns <- sample(ns, size = n.seed - length(gs), 
                      prob = net$combined_score[as.character(net$protein1) %in% 
                        gs])
                  gs <- unique(c(gs, ns))
                  qiter <- qiter + 1
                }
                rm[gs] <- k
                sampled[gs] <- TRUE
                if (n.seed <= 2) 
                  sampled[as.character(net$protein2[as.character(net$protein1) %in% 
                    gs])] <- TRUE
            }
            rm[attr(ratios, "rnames")]
        })
        rm <- do.call(cbind, tmp)
    }
    else if (substr(seed.method, 1, 7) == "netcor=") {
        if (!exists("ratios")) 
            stop("netcor seed method but no ratios")
        if (!exists("networks") || length(networks) <= 0) 
            stop("netcor seed method but no networks")
        seed.method <- strsplit(seed.method, "=")[[1]][2]
        net.name <- strsplit(seed.method, ":")[[1]][1]
        net <- networks[[net.name]]
        n.seed <- as.integer(strsplit(seed.method, ":")[[1]][2])
        rats <- get.cluster.matrix()
        cors <- cor(t(rats), use = "pairwise")
        tmp.mat <- matrix(0, nrow = nrow(cors), ncol = ncol(cors))
        dimnames(tmp.mat) <- dimnames(cors)
        tmp.lookup <- 1:attr(ratios, "nrow")
        names(tmp.lookup) <- attr(ratios, "rnames")
        net <- net[as.character(net$protein1) %in% attr(ratios, 
            "rnames") & as.character(net$protein2) %in% attr(ratios, 
            "rnames"), ]
        tmp.mat[cbind(tmp.lookup[as.character(net$protein1)], 
            tmp.lookup[as.character(net$protein2)])] <- net$combined_score/1000
        cors <- cors + tmp.mat
        rm(tmp.mat)
        rm <- rep(0, attr(ratios, "nrow"))
        names(rm) <- attr(ratios, "rnames")
        sampled <- rep(FALSE, attr(ratios, "nrow"))
        names(sampled) <- attr(ratios, "rnames")
        mc <- get.parallel(n.clust.per.row)
        tmp <- mc$apply(1:n.clust.per.row, function(i) {
            for (k in 1:k.clust) {
                if (sum(!sampled) < n.seed) 
                  sampled[sample(1:length(sampled))] <- FALSE
                rnames <- attr(ratios, "rnames")[!sampled]
                g <- sample(rnames, 1)
                g <- rnames[order(cors[g, !sampled], decreasing = T)[1:n.seed]]
                rm[g] <- k
                if (length(g) == 1) {
                  tmp <- rnames[order(cors[g, !sampled], decreasing = T)[1:10]]
                  sampled[tmp] <- TRUE
                }
                sampled[g] <- TRUE
            }
            rm[attr(ratios, "rnames")]
        })
        rm <- do.call(cbind, tmp)
    }
    if (is.vector(rm)) 
        rm <- t(rm)
    if (nrow(rm) == 1) 
        rm <- t(rm)
    if (col.method == "rnd") {
        cm <- t(sapply(1:attr(ratios, "ncol"), function(i) sample(1:k.clust, 
            n.clust.per.col[1], replace = n.clust.per.col[1] > 
                k.clust)))
    }
    else if (col.method == "best") {
        if (!exists("ratios")) 
            stop("best col seed method but no ratios")
        all.rats <- get.cluster.matrix()
        attr(all.rats, "all.colVars") <- apply(all.rats, 2, var, 
            use = "pair", na.rm = T)
        col.scores <- -sapply(1:k.clust, function(k) if (sum(rm == 
            k, na.rm = T) <= 0) 
            rep(NA, attr(ratios, "ncol"))
        else get.col.scores(k = rownames(which(rm == k, arr = T)), 
            ratios = all.rats, method = "orig"))
        cm <- t(apply(col.scores, 1, function(i) order(i, decreasing = T)[1:n.clust.per.col[1]]))
    }
    rownames(rm) <- attr(ratios, "rnames")
    rownames(cm) <- attr(ratios, "cnames")
    list(rm = rm, cm = cm)
}
set.param <-
function (name, val, env = cmonkey.params, override = F, quiet = F) 
{
    if (!exists(name, envir = env) || override) {
        if (!quiet) 
            try({
                cat(name, "-> ")
                str(val, digits.d = 9, no.list = T)
            })
        assign(name, val, envir = env)
    }
    else {
        val <- get(name, envir = env)
        if (!quiet) 
            try({
                cat(name, "= ")
                str(val, digits.d = 9, no.list = T)
            })
        assign(name, val, envir = env)
    }
    assign(name, val, envir = parent.frame())
}
spacer.one.cluster <-
function (k, seq.type = "upstream spacer", hits.to.all = F, verbose = F, 
    unlink = T, score.cutoff = -3, ...) 
{
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    bg.fname <- paste(progs.dir, "/SPACER/data/genomes/", organism, 
        ".upstream.raw", sep = "")
    if (!file.exists(bg.fname)) 
        dir.create(sprintf("%s/SPACER/data/genomes", progs.dir))
    writeLines(genome.info$all.upstream.seqs[[seq.type]], con = bg.fname)
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    filter <- FALSE
    if ("filter" %in% names(list(...)) && list(...)$filter && 
        !grepl("spacer", seq.type)[1]) 
        filter <- TRUE
    seqs <- get.sequences(rows, seq.type = seq.type, filter = filter, 
        ...)
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    bg.file <- paste(organism, ".upstream.raw", sep = "")
    tempfile <- my.tempfile(sprintf("spacer_%d_%d_", k, iter), 
        suf = ".fasta")
    temp.out <- my.tempfile(sprintf("spacer_out_%d_%d_", k, iter), 
        suf = ".txt")
    cat(paste(">", names(seqs), "\n", seqs, sep = ""), file = tempfile, 
        sep = "\n")
    cwd <- setwd(sprintf("%s/SPACER", progs.dir))
    on.exit(setwd(cwd))
    if (!grepl("prism", seq.type)[1]) {
        cmd <- sprintf(spacer.cmd[1], bg.file, temp.out, tempfile)
    }
    else {
        cmd <- sprintf(gsub("SPACER.jar", "PRISM.jar", spacer.cmd[1]), 
            bg.file, temp.out, tempfile)
    }
    if (verbose) 
        print(cmd)
    out <- system(cmd, intern = T, ignore.stderr = !verbose)
    if (!file.exists(temp.out)) 
        return(list(k = k))
    spacer.out <- readLines(temp.out)
    if (unlink) 
        unlink(temp.out)
    out <- strsplit(spacer.out[spacer.out != ""], "[\\t\\,]", 
        perl = T)
    all.fasta <- my.tempfile(sprintf("spacer_all_%d_%d_", k, 
        iter), suf = ".fasta")
    bg.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    if (hits.to.all && !is.na(bg.seqs) && length(out) > 1) 
        cat(paste(">", names(bg.seqs), "\n", bg.seqs, sep = ""), 
            file = all.fasta, sep = "\n")
    out2 <- out3 <- list()
    for (i in 1:length(out)) {
        motif <- out[[i]][2]
        motif.score <- as.numeric(out[[i]][1])
        if (motif.score < score.cutoff) 
            break
        prob <- 2^(-motif.score)
        temp.out <- my.tempfile(sprintf("spacer_out_%d_%d_", 
            k, iter), suf = ".txt")
        if (!grepl("prism", seq.type)[1]) {
            cmd <- sprintf(spacer.cmd[2], motif, temp.out, tempfile)
        }
        else {
            cmd <- sprintf(gsub("SPACER.jar", "PRISM.jar", spacer.cmd[2]), 
                motif, temp.out, tempfile)
        }
        if (verbose) 
            print(cmd)
        tmp <- system(cmd, intern = T, ignore.stderr = !verbose)
        tmp <- readLines(temp.out)
        spacer.out <- c(spacer.out, tmp)
        if (unlink && file.exists(temp.out)) 
            unlink(temp.out)
        tmp <- do.call(rbind, strsplit(tmp[tmp != ""], ","))
        out2[[length(out2) + 1]] <- tmp
        if (hits.to.all) {
            tmp.out2 <- my.tempfile(sprintf("spacer_out_%d_%d_", 
                k, iter))
            if (!file.exists(all.fasta)) {
                tmp <- bg.seqs[!bg.seqs %in% seqs]
                cat(paste(">", names(bg.seqs), "\n", bg.seqs, 
                  sep = ""), file = all.fasta, sep = "\n")
            }
            if (!grepl("prism", seq.type)[1]) {
                cmd <- sprintf(spacer.cmd[2], motif, tmp.out2, 
                  all.fasta)
            }
            else {
                cmd <- sprintf(gsub("SPACER.jar", "PRISM.jar", 
                  spacer.cmd[2]), motif, temp.out2, all.fasta)
            }
            if (verbose) 
                print(cmd)
            tmp <- system(cmd, intern = T, ignore.stderr = !verbose)
            tmp <- readLines(tmp.out2)
            tmp <- do.call(rbind, strsplit(tmp[tmp != ""], ","))
            tmp[, 1] <- names(bg.seqs)[as.integer(as.character(tmp[, 
                1]))]
            out3[[length(out3) + 1]] <- tmp
            if (unlink && file.exists(tmp.out2)) 
                unlink(tmp.out2)
        }
    }
    if (unlink && file.exists(tempfile)) 
        file.remove(tempfile)
    if (unlink && file.exists(all.fasta)) 
        file.remove(all.fasta)
    setwd(cwd)
    out <- list(motifs = do.call(rbind, out), hits = out2)
    attr(out, "spacer.out") <- spacer.out
    if (hits.to.all) 
        out$all.hits = out3
    meme.let <- c("A", "C", "G", "T")
    out$motifs <- data.frame(score = as.numeric(out$motifs[, 
        1]), motif = out$motifs[, 2], flag = out$motifs[, 3])
    out$motifs <- subset(out$motifs, score >= score.cutoff)
    if (length(out$hits) <= 0 || nrow(out$motifs) <= 0) 
        return(list(k = k, spacer.out = out))
    for (i in 1:length(out$hits)) {
        hits <- data.frame(seq = as.integer(out$hits[[i]][, 1]), 
            posn = as.integer(out$hits[[i]][, 2]), strand = out$hits[[i]][, 
                3], site = out$hits[[i]][, 4])
        sites <- toupper(do.call(rbind, strsplit(as.character(hits$site), 
            "")))
        if (any(!sites %in% meme.let)) 
            sites[!sites %in% meme.let] <- sample(meme.let, sum(!sites %in% 
                meme.let))
        pssm <- matrix(0, nrow = ncol(sites), ncol = 4)
        rownames(pssm) <- as.character(1:nrow(pssm))
        colnames(pssm) <- meme.let
        for (j in 1:nrow(sites)) pssm[cbind(1:nrow(pssm), sites[j, 
            ])] <- pssm[cbind(1:nrow(pssm), sites[j, ])] + 1
        out$hits[[i]] <- list(sites = hits, pssm = pssm)
    }
    m.in <- character()
    for (i in 1:length(out$hits)) {
        if (out$motifs$score[i] < score.cutoff) 
            next
        m.in <- c(m.in, pssm.motif.lines(out$hits[[i]]$pssm, 
            id = sprintf("spacer_%d", i), header = (i == 1)))
    }
    all.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    mast.out <- runMast(m.in, mast.cmd[seq.type], names(all.seqs), 
        all.seqs, bg.list = genome.info$bg.list[[seq.type]], 
        unlink = T, verbose = verbose)
    pv.ev <- get.pv.ev.single(mast.out, rows)
    meme.out <- list()
    for (ii in 1:length(out$hits)) {
        wo <- out$hits[[ii]]
        pssm <- wo$pssm
        pssm <- pssm + max(pssm, na.rm = T)/100
        for (i in 1:nrow(pssm)) pssm[i, ] <- pssm[i, ]/sum(pssm[i, 
            ], na.rm = T)
        posns <- data.frame(gene = names(seqs)[wo$sites$seq], 
            strand = wo$sites$strand, start = as.integer(wo$sites$posn), 
            p.value = NA, site = toupper(wo$sites$site))
        meme.out[[ii]] <- list(width = nrow(wo$pssm), sites = nrow(wo$sites), 
            llr = out$motifs$score[ii], e.value = 2^(-out$motifs$score[ii]), 
            pssm = pssm, posns = posns)
    }
    attr(meme.out, "is.pal") <- FALSE
    invisible(list(k = k, spacer.out = out, meme.out = meme.out, 
        pv.ev = pv.ev))
}
system.time.limit <-
function (cmd, tlimit = 600) 
{
    out <- readLines(pipe(cmd, open = "rt"))
    out
}
test.fit.cluster <-
function (k, verbose = F, plot = F, ...) 
{
    scores <- get.all.scores(k, verbose, ...)
    x <- scores$r
    row.memb <- attr(ratios, "rnames") %in% get.rows(k)
    names(row.memb) <- attr(ratios, "rnames")
    y <- as.factor(row.memb)
    xx <- x[, apply(x[y == TRUE, ], 2, sum) != 0]
    if (plot) {
        len <- 25
        xp <- expand.grid(as.data.frame(apply(x, 2, function(i) seq(min(i), 
            max(i), length = len))))
        xxp <- expand.grid(as.data.frame(apply(xx, 2, function(i) seq(min(i), 
            max(i), length = len))))
    }
    out0 <- glm(y ~ . - 1, data = as.data.frame(x), family = "binomial", 
        ...)
    if (verbose) 
        print(summary(out0))
    prob0 <- predict(out0, type = "response")
    names(prob0) <- row.names(row.scores)
    if (plot) 
        prob0p <- predict(out0, newdata = as.data.frame(xp), 
            type = "response")
    require(varSelRF)
    out1 <- varSelRF(x, y, ntree = 5000, ntreeIterat = 2000, 
        vars.drop.frac = 0.5, whole.range = F, keep.forest = T, 
        ...)
    if (verbose) 
        print(out1)
    prob1 <- predict(out1$rf.model, type = "prob", newdata = subset(x, 
        select = out1$selected.vars))[, 2]
    if (plot) 
        prob1p <- predict(out1$rf.model, type = "prob", newdata = subset(xp, 
            select = out1$selected.vars))[, 2]
    require(brglm)
    out2 <- brglm(y ~ . - 1, data = as.data.frame(x), ...)
    if (verbose) 
        print(summary(out2))
    prob2 <- predict(out2, type = "response")
    names(prob2) <- row.names(row.scores)
    if (plot) 
        prob2p <- predict(out2, newdata = as.data.frame(xp), 
            type = "response")
    require(nnet)
    out3 <- nnet(y ~ . - 1, data = as.data.frame(x), skip = T, 
        softmax = F, size = 3, decay = 0.1, maxit = 1000, tra = F, 
        ...)
    prob3 <- predict(out3)[, 1]
    if (plot) 
        prob3p <- predict(out3, newdata = as.data.frame(xp))
    require(MASS)
    prob4 <- prob5 <- rep(NA, length(prob1))
    out4 <- try(qda(y ~ . - 1, data = as.data.frame(xx), method = "mle", 
        ...))
    if (!"try-error" %in% class(out4)) {
        prob4 <- predict(out4)$posterior[, 2]
        if (plot) 
            prob4p <- predict(out4, newdata = as.data.frame(xp))$posterior[, 
                2]
    }
    out5 <- try(lda(y ~ . - 1, data = as.data.frame(xx), method = "mle", 
        ...))
    if (!"try-error" %in% class(out5)) {
        prob5 <- predict(out5)$posterior[, 2]
        if (plot) 
            prob5p <- predict(out5, newdata = as.data.frame(xp))$posterior[, 
                2]
    }
    require(class)
    out6 <- knn.cv(x, as.integer(y), k = attr(ratios, "nrow")/k.clust * 
        n.clust.per.row, prob = T, use.all = T, ...)
    prob6 <- as.numeric(attr(out6, "prob"))
    if (plot) {
        prob6p <- knn(x, xp, as.integer(y), k = attr(ratios, 
            "nrow")/k.clust * n.clust.per.row, prob = T, use.all = T, 
            ...)
        prob6p <- as.numeric(attr(prob6p, "prob"))
    }
    require(e1071)
    out7 <- svm(y ~ . - 1, data = as.data.frame(cbind(x, y = as.integer(y) - 
        1)), probability = T, scale = F, ...)
    prob7 <- predict(out7, as.data.frame(x), prob = T)
    if (plot) 
        prob7p <- predict(out7, as.data.frame(xp), prob = T)
    require(mda)
    out8 <- fda(y ~ . - 1, data = as.data.frame(cbind(x, y = as.integer(y) - 
        1)), ...)
    prob8 <- predict(out8, type = "posterior")[, 2]
    if (plot) 
        prob8p <- predict(out8, as.data.frame(xp), type = "posterior")[, 
            2]
    out9 <- mda(y ~ . - 1, data = as.data.frame(cbind(x, y = as.integer(y) - 
        1)), ...)
    prob9 <- predict(out9, newdata = as.data.frame(x), type = "posterior")[, 
        2]
    if (plot) 
        prob9p <- predict(out9, as.data.frame(xp), type = "posterior")[, 
            2]
    probs <- cbind(prob0, prob1, prob2, prob3, prob4, prob5, 
        prob6, prob7, prob8, prob9)
    if (plot) {
        xp2 <- apply(x, 2, function(i) seq(min(i), max(i), length = len))
        probsp <- cbind(prob0p, prob1p, prob2p, prob3p, prob4p, 
            prob5p, prob6p, prob7p, prob8p, prob9p)
        probsp.tot <- array(apply(probsp, 1, mean, na.rm = T), 
            dim = rep(len, ncol(xp)))
        comb <- t(combn(1:ncol(x), 2))
        par(mfrow = c(ceiling(sqrt(nrow(comb))), floor(sqrt(nrow(comb)))))
        for (i in 1:nrow(comb)) {
            i1 <- comb[i, 1]
            i2 <- comb[i, 2]
            plot(x[, i1], x[, i2], xlab = colnames(x)[i1], ylab = colnames(x)[i2], 
                col = as.integer(y), pch = 20, cex = 0.5)
            tmp <- aperm(probsp.tot, c(i1, i2, (1:ncol(x))[!1:ncol(x) %in% 
                comb[i, ]]))
            contour(xp2[, i1], xp2[, i2], tmp[, , 1, 1], add = T, 
                labex = 0)
        }
    }
    invisible(probs)
}
un.ffify.env <-
function (env) 
{
    for (i in c("row.scores", "mot.scores", "net.scores")) if (exists(i, 
        envir = env)) 
        env[[i]] <- env[[i]][, ]
    for (i in names(env$ratios)) if (!is.null(env$ratios[[i]])) 
        env$ratios[[i]] <- env$ratios[[i]][, ]
    for (i in names(env$meme.scores)) if (!is.null(env$meme.scores[[i]])) 
        env$meme.scores[[i]] <- as.list(env$meme.scores[[i]])
    for (i in c("clusterStack", "genome.info", "networks")) if (exists(i, 
        envir = env)) 
        env[[i]] <- as.list(env[[i]])
    invisible(env)
}
update.all.clusters <-
function (env, dont.update = F, ...) 
{
    mc <- get.parallel(k.clust)
    all.scores <- mc$apply(1:k.clust, get.all.scores, return.scores = F, 
        densities = T, verbose = F, members = T, force.motif = F, 
        ...)
    if (!dont.update) {
        for (k in 1:length(all.scores)) {
            if (!is.null(all.scores[[k]]$ms)) {
                for (i in names(all.scores[[k]]$ms)) if (!is.null(all.scores[[k]]$ms[[i]])) 
                  env$meme.scores[[i]][[k]] <- all.scores[[k]]$ms[[i]]
            }
            env$clusterStack[[k]]$rows <- all.scores[[k]]$members$r
            env$clusterStack[[k]]$nrows <- length(all.scores[[k]]$members$r)
            env$clusterStack[[k]]$cols <- all.scores[[k]]$members$c
            env$clusterStack[[k]]$ncols <- length(all.scores[[k]]$members$c)
        }
        env$clusterStack <- env$get.clusterStack(ks = 1:k.clust, 
            force = T)
    }
    env
}
update.cmonkey.env <-
function (object, ...) 
{
    if (file.exists("cmonkey-funcs.R")) {
        tmp.e <- new.env()
        sys.source("cmonkey.R", envir = tmp.e)
    }
    else {
        tmp.e <- environment(cMonkey:::cmonkey)
    }
    for (i in ls(tmp.e)) {
        if (i %in% c("DATE", "VERSION")) 
            next
        f <- try(get(i, envir = tmp.e))
        f2 <- try(get(paste("super", i, sep = "."), envir = object), 
            silent = T)
        if (class(f) == "function") {
            environment(f) <- object
            if (class(f2) != "function") 
                assign(i, f)
            else assign(paste("super", i, sep = "."), f)
        }
    }
    rm(f, f2, tmp.e, i)
    for (i in ls()) {
        if (i %in% c("i", "object")) 
            next
        f <- get(i)
        if (is.function(f)) 
            assign(i, f, object)
    }
    if (FALSE && (env$big.memory == TRUE || env$big.memory > 
        0)) {
        for (i in c("row.scores", "mot.scores", "net.scores", 
            "col.scores")) if (!is.null(env[[i]]) && require(bigmemory) && 
            is.big.matrix(env[[i]])) 
            attach.big.matrix(get(i, env))
        for (i in 1:length(env$ratios)) {
            if (!is.null(env$ratios[[i]]) && require(bigmemory) && 
                is.big.matrix(env$ratios[[i]])) 
                attach.big.matrix(env$ratios[[i]])
        }
    }
}
viewPssm <-
function (pssm, e.val = NA, mot.ind = NA, use.char = T, main.title = NA, 
    no.par = F, ...) 
{
    if (is.null(pssm)) 
        return()
    getEntropy <- function(pssm) {
        pssm[pssm == 0] <- 1e-05
        entropy <- apply(pssm, 1, function(i) -sum(i * log2(i)))
        return(entropy)
    }
    char.coords = list(T = list(x = c(0.45, 0.55, 0.55, 1, 1, 
        0, 0, 0.45), y = c(0, 0, 0.9, 0.9, 1, 1, 0.9, 0.9), color = 2), 
        A = list(x = c(0, 0.1, 0.28, 0.72, 0.68, 0.32, 0.5, 0.9, 
            1, 0.55, 0.45, 0), y = c(0, 0, 0.4, 0.4, 0.5, 0.5, 
            0.9, 0, 0, 1, 1, 0), color = 3), C = list(x = c(1, 
            1, 0.85, 0.55, 0.45, 0.15, 0, 0, 0.15, 0.45, 0.55, 
            0.85, 1, 1, 0.9, 0.9, 0.8, 0.55, 0.45, 0.2, 0.1, 
            0.1, 0.2, 0.45, 0.55, 0.8, 0.9, 0.9), y = c(0.6, 
            0.7, 0.9, 1, 1, 0.9, 0.65, 0.35, 0.1, 0, 0, 0.1, 
            0.35, 0.4, 0.4, 0.35, 0.2, 0.1, 0.1, 0.2, 0.42, 0.58, 
            0.8, 0.9, 0.9, 0.8, 0.65, 0.6), color = 4), G = list(x = c(1, 
            1, 0.85, 0.55, 0.45, 0.15, 0, 0, 0.15, 0.45, 0.55, 
            0.85, 1, 1, 0.7, 0.7, 0.9, 0.8, 0.55, 0.45, 0.2, 
            0.1, 0.1, 0.2, 0.45, 0.55, 0.8, 0.9, 0.9), y = c(0.6, 
            0.7, 0.9, 1, 1, 0.9, 0.65, 0.35, 0.1, 0, 0, 0.1, 
            0.35, 0.5, 0.5, 0.4, 0.4, 0.2, 0.1, 0.1, 0.2, 0.42, 
            0.58, 0.8, 0.9, 0.9, 0.8, 0.65, 0.6), color = "orange"))
    draw.char <- function(char = col.let, rect = c(0, 0, 1, 1)) {
        if (rect[4] <= 1e-05) 
            return()
        x <- char.coords[[char]]$x * rect[3] + rect[1]
        y <- char.coords[[char]]$y * rect[4] + rect[2]
        color <- char.coords[[char]]$color
        polygon(x, y, col = color, border = color)
    }
    win.size <- nrow(pssm)
    if (!no.par) 
        par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 0) * 0.75)
    if (any(pssm <= 0)) 
        pssm <- pssm + 1e-10
    if (any(pssm > 1)) 
        pssm <- t(apply(pssm, 1, function(i) i/(sum(i) + 1e-10)))
    entr <- getEntropy(pssm)
    scale.e <- (2 - entr)/2
    x.range <- c(0.5, win.size + 0.5)
    y.range <- c(0, 1)
    plot(x.range, y.range, type = "n", tck = 0.01, cex.lab = 0.2, 
        cex.sub = 0.2, cex.axis = 0.2, axes = F)
    if (!is.na(main.title[1])) {
        if (!is.na(mot.ind)) 
            title(main.title, col.main = mot.ind + 1, xpd = NA, 
                ...)
        else title(main.title, xpd = NA, ...)
    }
    else if (!is.na(mot.ind) || !is.na(e.val)) {
        if (!is.na(mot.ind)) {
            if (!is.na(e.val)) 
                tmp.tit <- sprintf("PSSM #%d; E=%.3g", mot.ind, 
                  e.val)
            else tmp.tit <- tmp.tit <- sprintf("PSSM #%d", mot.ind)
            title(tmp.tit, col.main = mot.ind + 1, xpd = NA, 
                ...)
        }
        else title(tmp.tit, xpd = NA, ...)
    }
    pssm.sc <- scale.e * pssm
    for (j in 1:win.size) {
        inds <- sort(pssm.sc[j, ], index = T)$ix
        for (i in 1:4) {
            ind <- inds[i]
            if (i == 1) {
                if (!use.char) {
                  rect((j - 0.5), 0, (j + 0.5), pssm.sc[j, ind], 
                    col = colMap[ind])
                  if (pssm[j, ind] > 0.05) 
                    text(j, 0 + pssm.sc[j, ind]/2, colLet[ind])
                }
                else {
                  draw.char(col.let[ind], c((j - 0.4), 0, 0.9, 
                    pssm.sc[j, ind] - 0.01))
                }
                prev.h <- pssm.sc[j, ind]
            }
            else {
                if (!use.char) {
                  rect((j - 0.5), prev.h, (j + 0.5), (pssm.sc[j, 
                    ind] + prev.h), col = colMap[ind])
                  if (pssm.sc[j, ind] > 0.05) {
                    if (i == 2) 
                      text(j, prev.h + 0.5 * pssm.sc[j, ind], 
                        colLet[ind], col = 8)
                    else text(j, prev.h + 0.5 * pssm.sc[j, ind], 
                      colLet[ind])
                  }
                }
                else {
                  draw.char(col.let[ind], c((j - 0.4), prev.h, 
                    0.9, pssm.sc[j, ind] - 0.01))
                }
                prev.h <- prev.h + pssm.sc[j, ind]
            }
        }
        if (win.size < 10) 
            text(1:win.size, rep(-0.01, win.size), as.character(1:win.size), 
                cex = 0.7, adj = c(0.5, 1), xpd = NA)
        else if (win.size < 20) 
            text(seq(1, win.size, 2), rep(-0.01, win.size), as.character(seq(1, 
                win.size, 2)), cex = 0.7, adj = c(0.5, 1), xpd = NA)
        else text(seq(1, win.size, 5), rep(-0.01, win.size), 
            as.character(seq(1, win.size, 5)), cex = 0.7, adj = c(0.5, 
                1), xpd = NA)
    }
}
weeder.one.cluster <-
function (k, seq.type = "upstream weeder", n.motifs = 4, verbose = F, 
    unlink = T, weeder.size = "medium", ...) 
{
    ntides <- c("T", "G", "A", "C")
    for (w in c(6, 8)) {
        if (!file.exists(sprintf("%s/FreqFiles", progs.dir))) 
            dir.create(sprintf("%s/FreqFiles", progs.dir))
        if (!file.exists(paste(sprintf("%s/FreqFiles/", progs.dir), 
            toupper(organism), ".", w, ".", "freq", sep = ""))) {
            seqs <- unique(genome.info$all.upstream.seqs[[seq.type]])
            all.substrings <- as.vector(sapply(1:(max(nchar(seqs)) - 
                w + 1), function(i) substr(seqs, i, i + w - 1)))
            all.substrings <- all.substrings[!is.na(all.substrings) & 
                all.substrings != "" & nchar(all.substrings) == 
                w]
            all.substrings <- all.substrings[!grepl("[^GATC]", 
                all.substrings)]
            hist.substrings <- table(as.factor(all.substrings))
            all.combos <- all.dna.seqs(w, ntides)
            all.combos <- all.combos[!all.combos %in% names(hist.substrings)]
            tmp <- rep(0, length(all.combos))
            names(tmp) <- all.combos
            hist.substrings <- c(hist.substrings, tmp) + 1
            hist.substrings <- hist.substrings[sort(names(hist.substrings))]
            write.table(hist.substrings, quote = F, sep = " ", 
                col.names = F, file = paste(sprintf("%s/FreqFiles/", 
                  progs.dir), toupper(organism), ".", w, ".", 
                  "freq", sep = ""))
        }
    }
    if (is.numeric(k)) 
        rows <- get.rows(k)
    else rows <- k
    seqs <- get.sequences(rows, seq.type = seq.type, ...)
    min.seqs <- cluster.rows.allowed[1]
    max.seqs <- cluster.rows.allowed[2]
    if (is.null(seqs) || length(seqs) < min.seqs) 
        return(list(k = k))
    if (length(seqs) < min.seqs || length(seqs) > max.seqs) 
        return(list(k = k))
    cat(k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", 
        length(seqs), "\n")
    fst.file <- my.tempfile(sprintf("weeder_%d_%d_", k, iter), 
        suf = ".fst")
    cat(paste(">", names(seqs), "\n", seqs, sep = ""), file = fst.file, 
        sep = "\n")
    file.remove(paste(fst.file, c("fst", "html", "mix", "wee"), 
        sep = "."))
    cwd <- setwd(progs.dir)
    on.exit(setwd(cwd))
    cmd <- sprintf(weeder.cmd, fst.file, toupper(organism), weeder.size, 
        n.motifs * 5)
    if (verbose) 
        print(cmd)
    out <- system(cmd, intern = T, ignore.stderr = !verbose)
    setwd(cwd)
    if (!file.exists(paste(fst.file, "mix", sep = ".")) && !file.exists(paste(fst.file, 
        "wee", sep = "."))) 
        return(list(k = k))
    out <- c(readLines(paste(fst.file, "mix", sep = ".")), readLines(paste(fst.file, 
        "wee", sep = ".")))
    if (unlink) 
        file.remove(c(fst.file, paste(fst.file, c("fst", "html", 
            "mix", "wee"), sep = ".")))
    mot.scores <- as.data.frame(do.call(rbind, strsplit(grep("^\\d+\\) ", 
        out, perl = T, val = T), " "))[, 2:4])
    mot.scores[, 2] <- as.numeric(as.character(mot.scores[, 2]))
    mot.scores[, 3] <- as.numeric(as.character(mot.scores[, 3]))
    mot.scores[, 3][is.na(mot.scores[, 3])] <- 0
    is.highest.ranking.motif <- grep("Interesting motifs (highest-ranking) seem to be", 
        out, fixed = T)
    is.not.highest.ranking.motif <- grep("Interesting motifs (not highest-ranking) can also be", 
        out, fixed = T)
    starts <- grep("Best occurrences", out, fixed = T)
    weeder.out <- list()
    for (i in 1:length(starts)) {
        start <- starts[i]
        weeder.out[[i]] <- list()
        weeder.out[[i]]$motifs.redund <- strsplit(out[start - 
            2], "\\s\\-\\s", perl = T)[[1]]
        weeder.out[[i]]$motifs <- out[c(start - 5, start - 6)]
        weeder.out[[i]]$is.highest.ranking <- start > is.highest.ranking.motif && 
            start < is.not.highest.ranking.motif
        weeder.out[[i]]$score <- unique(subset(mot.scores, V1 %in% 
            weeder.out[[i]]$motifs)$V2)
        end <- start - 1 + min(which(out[(start + 1):length(out)] == 
            ""))
        lines <- strsplit(out[(start + 1):end], "\\s+", perl = T)
        lines <- do.call(rbind, lines)
        colnames(lines)[2:ncol(lines)] <- lines[1, 1:(ncol(lines) - 
            1)]
        lines <- as.data.frame(lines[-1, -1])
        lines$match <- gsub("(", "", gsub(")", "", lines$match, 
            fixed = T), fixed = T)
        weeder.out[[i]]$matches <- lines
        i2 <- end + min(grep("Frequency Matrix", out[end:length(out)])) + 
            1
        lines <- strsplit(out[i2:(i2 - 1 + min(which(out[(i2 + 
            1):length(out)] == "")))], "\t+", perl = T)
        lines <- do.call(rbind, lines)[, -1]
        counts1 <- do.call(rbind, strsplit(lines[, 1], "\\s+", 
            perl = T))[, -1]
        colnames(counts1) <- counts1[1, ]
        counts1 <- counts1[-1, ]
        counts1 <- apply(counts1, 2, as.integer)
        counts2 <- do.call(rbind, strsplit(lines[, 2], "\\s+", 
            perl = T))[, -1]
        colnames(counts2) <- counts2[1, ]
        counts2 <- counts2[-1, ]
        counts2 <- apply(counts2, 2, as.integer)
        weeder.out[[i]]$counts.all <- counts1
        weeder.out[[i]]$counts.best <- counts2
    }
    if (length(weeder.out) <= 0) {
        attr(weeder.out, "weeder.out") <- out
        return(list(k = k, weeder.out = weeder.out))
    }
    weeder.out <- weeder.out[order(sapply(weeder.out, function(i) nchar(i$motifs[1])))]
    n.mot <- length(weeder.out)
    if (n.mot > 1) {
        m.redund <- matrix(0, nrow = n.mot, ncol = n.mot)
        for (i in 1:(n.mot - 1)) for (j in (i + 1):n.mot) m.redund[i, 
            j] <- sum(weeder.out[[i]]$motifs %in% weeder.out[[j]]$motifs.redund)
        weeder.out <- weeder.out[order(apply(m.redund, 1, sum, 
            na.rm = T), decreasing = F)]
        m.redund <- m.redund * 0
        n.mot <- length(weeder.out)
        for (i in 1:(n.mot - 1)) for (j in (i + 1):n.mot) m.redund[i, 
            j] <- sum(weeder.out[[i]]$motifs %in% weeder.out[[j]]$motifs.redund)
        sum.redund <- apply(m.redund, 1, sum, na.rm = T)
        m.length <- sapply(weeder.out, function(i) nchar(i$motifs[1]))
        weeder.out <- weeder.out[order(sapply(weeder.out, "[[", 
            "is.highest.ranking"), m.length, sum.redund, decreasing = T)]
        m.redund <- m.redund * 0
        n.mot <- length(weeder.out)
        for (i in 1:(n.mot - 1)) for (j in (i + 1):n.mot) m.redund[i, 
            j] <- sum(weeder.out[[i]]$motifs %in% weeder.out[[j]]$motifs.redund)
        attr(weeder.out, "is.redund") <- m.redund
        weeder.out <- weeder.out[1:n.motifs]
        weeder.out <- weeder.out[!sapply(weeder.out, is.null)]
    }
    attr(weeder.out, "weeder.out") <- out
    m.in <- character()
    for (i in 1:length(weeder.out)) m.in <- c(m.in, pssm.motif.lines(weeder.out[[i]]$counts.all, 
        id = sprintf("weeder_%d", i), header = (i == 1)))
    all.seqs <- genome.info$all.upstream.seqs[[seq.type]]
    mast.out <- runMast(m.in, mast.cmd[seq.type], names(all.seqs), 
        all.seqs, bg.list = genome.info$bg.list[[seq.type]], 
        unlink = T, verbose = verbose)
    pv.ev <- get.pv.ev.single(mast.out, rows)
    meme.out <- list()
    for (ii in 1:length(weeder.out)) {
        wo <- weeder.out[[ii]]
        pssm <- wo$counts.all
        pssm <- pssm + max(pssm, na.rm = T)/100
        for (i in 1:nrow(pssm)) pssm[i, ] <- pssm[i, ]/sum(pssm[i, 
            ], na.rm = T)
        posns <- data.frame(gene = names(seqs)[wo$matches$Seq], 
            strand = wo$matches$St, start = as.integer(wo$matches$pos), 
            p.value = (100 - as.numeric(wo$matches$match) + 0.001)/100, 
            site = gsub("[\\[\\]]", "", wo$matches$oligo, perl = T))
        meme.out[[ii]] <- list(width = nrow(wo$counts.all), sites = nrow(wo$matches), 
            llr = wo$score, e.value = wo$score, pssm = pssm, 
            posns = posns)
    }
    attr(meme.out, "is.pal") <- FALSE
    invisible(list(k = k, weeder.out = weeder.out, meme.out = meme.out, 
        pv.ev = pv.ev))
}
write.bicluster.network <-
function (out.dir = NULL, ks = 1:k.clust, seq.type = names(mot.weights)[1], 
    tomtom = T, tt.filter = function(tt) subset(tt, overlap >= 
        4 & q.value <= 0.05), m.filter = function(m) subset(m, 
        e.value <= Inf), gene.url = function(g) sprintf("http://microbesonline.org/cgi-bin/keywordSearch.cgi?taxId=%d&keyword=%s", 
        taxon.id, g), image.urls = T, ...) 
{
    if (is.null(out.dir)) {
        out.dir <- paste(cmonkey.filename, "network", sep = "/")
        if (iter != n.iter) 
            out.dir <- sprintf("%s_%04d/network", cmonkey.filename, 
                iter)
    }
    if (!file.exists(out.dir)) 
        dir.create(out.dir, recursive = T, showWarnings = F)
    cat("Outputing to", out.dir, "\n")
    r.sif <- do.call(rbind, lapply(ks, function(k) data.frame(get.rows(k), 
        sprintf("bicluster_%04d", k))))
    r.sif <- data.frame(r.sif[, 1], "gene_member", r.sif[, 2])
    colnames(r.sif) <- c("node1", "int", "node2")
    ms <- meme.scores[[seq.type]]
    m.sif <- do.call(rbind, lapply(ks, function(k) if (length(ms[[k]]$meme.out) <= 
        0) 
        NULL
    else data.frame(sprintf("motif_%04d_%d", k, 1:length(ms[[k]]$meme.out)), 
        sprintf("bicluster_%04d", k), sapply(ms[[k]]$meme.out, 
            "[[", "width"), sapply(ms[[k]]$meme.out, "[[", "sites"), 
        sapply(ms[[k]]$meme.out, "[[", "llr"), sapply(ms[[k]]$meme.out, 
            "[[", "e.value"))))
    m.sif <- data.frame(m.sif[, 1], "motif", m.sif[, 2:ncol(m.sif)])
    colnames(m.sif) <- c("node1", "int", "node2", "width", "nsites", 
        "llr", "e.value")
    if (!is.null(m.filter)) 
        m.sif <- m.filter(m.sif)
    if ("tout" %in% names(list(...))) 
        tout <- list(...)$tout
    if (!exists("tout")) 
        tout <- motif.similarities.tomtom(ks, ks, ...)
    tout <- subset(tout, !is.na(biclust1))
    tt.sif <- data.frame(sprintf("motif_%04d_%d", as.integer(tout[, 
        1]), as.integer(tout[, 2])), sprintf("motif_%04d_%d", 
        as.integer(tout[, 5]), as.integer(tout[, 6])), tout[, 
        9], tout[, 10], tout[, 11], tout[, 12])
    tt.sif <- data.frame(tt.sif[, 1], "motif_sim", tt.sif[, 2:ncol(tt.sif)])
    colnames(tt.sif) <- c("node1", "int", "node2", "offset", 
        "p.value", "q.value", "overlap")
    if (!is.null(tt.filter)) 
        tt.sif <- tt.filter(tt.sif)
    out.sif <- rbind(r.sif[, 1:3], m.sif[, 1:3], tt.sif[, 1:3])
    m.eda <- data.frame(edge = paste(m.sif[, 1], " (", m.sif[, 
        2], ") ", m.sif[, 3], sep = ""), m.sif[, 4:ncol(m.sif)])
    tt.eda <- data.frame(edge = paste(tt.sif[, 1], " (", tt.sif[, 
        2], ") ", tt.sif[, 3], sep = ""), tt.sif[, 4:ncol(tt.sif)])
    node.type <- as.data.frame(rbind(as.matrix(data.frame(attr(ratios, 
        "rnames"), "gene")), as.matrix(data.frame(attr(ratios, 
        "cnames"), "condition")), as.matrix(data.frame(sprintf("bicluster_%04d", 
        ks), "bicluster")), as.matrix(data.frame(as.character(m.sif$node1), 
        "motif"))))
    colnames(node.type) <- c("node", "type")
    syn.names <- do.call(rbind, lapply(attr(ratios, "rnames"), 
        function(g) data.frame(g, paste(get.synonyms(g)[[g]], 
            collapse = "::"))))
    colnames(syn.names) <- c("gene", "synonyms")
    l.names <- do.call(rbind, lapply(attr(ratios, "rnames"), 
        function(g) data.frame(g, paste(get.long.names(g)[[g]], 
            collapse = "::"))))
    colnames(l.names) <- c("gene", "long.name")
    s.names <- do.call(rbind, lapply(attr(ratios, "rnames"), 
        function(g) data.frame(g, paste(get.long.names(g, short = T)[[g]], 
            collapse = "::"))))
    colnames(s.names) <- c("gene", "short.name")
    g.url <- NULL
    if (!is.null(gene.url)) {
        g.url <- do.call(rbind, lapply(attr(ratios, "rnames"), 
            function(g) data.frame(g, gene.url(g))))
        colnames(g.url) <- c("gene", "url")
    }
    mot.info <- do.call(rbind, lapply(ks, function(k) if (length(ms[[k]]$meme) <= 
        0) 
        NULL
    else data.frame(sprintf("motif_%04d_%d", k, 1:length(ms[[k]]$meme.out)), 
        sapply(ms[[k]]$meme.out, function(i) pssm.to.string(i$pssm)), 
        sapply(ms[[k]]$meme.out, "[[", "width"), sapply(ms[[k]]$meme.out, 
            "[[", "sites"), sapply(ms[[k]]$meme.out, "[[", "llr"), 
        sapply(ms[[k]]$meme.out, "[[", "e.value"), sprintf("file://%s/%s/htmls/cluster%04d_pssm%d.png", 
            getwd(), out.dir, k, 1:length(ms[[k]]$meme.out)))))
    colnames(mot.info) <- c("motif", "consensus", "width", "n.sites", 
        "llr", "e.value", "imgURL")
    clust.info <- lapply(c("resid", "p.clust", "e.val", "nrows", 
        "ncols"), function(i) sapply(ks, function(j) clusterStack[[j]][[i]]))
    names(clust.info) <- c("resid", "p.clust", "e.val", "nrows", 
        "ncols")
    clust.info[[names(unlist(sapply(clust.info, nrow)))]] <- t(clust.info[[names(unlist(sapply(clust.info, 
        nrow)))]])
    for (n in names(clust.info)) if (!is.null(ncol(clust.info[[n]])) && 
        is.null(colnames(clust.info[[n]]))) 
        colnames(clust.info[[n]]) <- paste(n, 1:ncol(clust.info[[n]]), 
            sep = ".")
    clust.info <- cbind(bicluster = sprintf("bicluster_%04d", 
        ks), do.call(cbind, clust.info), url = sprintf("file://%s/%s/htmls/cluster%04d.html", 
        getwd(), out.dir, ks), imgURL = sprintf("file://%s/%s/htmls/cluster%04d_profile.png", 
        getwd(), out.dir, ks))
    rownames(clust.info) <- NULL
    noa <- merge(node.type, syn.names, by.x = "node", by.y = "gene", 
        all = T)
    noa <- merge(noa, l.names, by.x = "node", by.y = "gene", 
        all = T)
    noa <- merge(noa, s.names, by.x = "node", by.y = "gene", 
        all = T)
    noa <- merge(noa, mot.info, by.x = "node", by.y = "motif", 
        all = T)
    noa <- merge(noa, clust.info, by.x = "node", by.y = "bicluster", 
        all = T)
    if (!is.null(g.url)) 
        noa <- merge(noa, g.url, by.x = "node", by.y = "gene", 
            all = T)
    noa$imgURL <- as.character(noa$imgURL.x)
    noa$imgURL[!is.na(noa$imgURL.y)] <- as.character(noa$imgURL.y[!is.na(noa$imgURL.y)])
    noa <- noa[, !colnames(noa) %in% c("imgURL.x", "imgURL.y")]
    write.table(out.sif, quote = F, sep = "\t", col.names = T, 
        row.names = F, file = paste(out.dir, "all.sif", sep = "/"))
    write.table(noa, quote = F, sep = "\t", col.names = T, row.names = F, 
        file = paste(out.dir, "all.noa", sep = "/"))
    write.table(m.eda, quote = F, sep = "\t", col.names = T, 
        row.names = F, file = paste(out.dir, "m.eda", sep = "/"))
    write.table(tt.eda, quote = F, sep = "\t", col.names = T, 
        row.names = F, file = paste(out.dir, "tt.eda", sep = "/"))
    cat("Wrote", nrow(noa), "nodes and", nrow(out.sif), "edges to", 
        out.dir, "\n")
    invisible(list(sif = out.sif, noa = noa, m.eda = m.eda, tt.eda = tt.eda))
}
write.project <-
function (ks = sapply(as.list(clusterStack), "[[", "k"), para.cores = 1, 
    out.dir = NULL, gaggle = T, seq.type = names(mot.weights)[1], 
    gzip = T, output = c("svg", "pdf", "png", "html", "main", 
        "rdata"), ...) 
{
    if (is.null(out.dir)) {
        out.dir <- cmonkey.filename
        if (iter != n.iter) 
            out.dir <- sprintf("%s_%04d", out.dir, iter)
    }
    cat("Outputing to", out.dir, "\n")
    if (!file.exists(out.dir)) 
        dir.create(out.dir, recursive = T, showWarnings = F)
    clusterStack <- clusterStack[ks]
    mc <- get.parallel(length(ks), para.cores = para.cores)
    has.pdftk <- length(system("which pdftk", intern = T)) > 
        0
    if (!file.exists(paste(out.dir, "/svgs", sep = ""))) 
        dir.create(paste(out.dir, "/svgs", sep = ""), showWarnings = F)
    if (!file.exists(paste(out.dir, "/pdfs", sep = ""))) 
        dir.create(paste(out.dir, "/pdfs", sep = ""), showWarnings = F)
    if (!file.exists(paste(out.dir, "/htmls", sep = ""))) 
        dir.create(paste(out.dir, "/htmls", sep = ""), showWarnings = F)
    if ("svg" %in% output) {
        require(RSVGTipsDevice)
        if (!file.exists(sprintf("%s/svgs/stats.svg", out.dir)) && 
            !file.exists(sprintf("%s/svgs/stats.svgz", out.dir))) {
            cat("STATS...\n")
            devSVGTips(sprintf("%s/svgs/stats.svg", out.dir), 
                toolTipMode = 2, title = "Biclustering statistics", 
                xmlHeader = T)
            par(family = "Arial")
            plotStats(new.dev = F)
            dev.off()
        }
        require(igraph)
        cat("SVGS: ")
        for (qqq in 1:3) {
            lapply(ks, function(k) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (file.exists(sprintf("%s/svgs/cluster%04d.svg", 
                  out.dir, k)) || file.exists(sprintf("%s/svgs/cluster%04d.svgz", 
                  out.dir, k))) 
                  return(NULL)
                devSVGTips(sprintf("%s/svgs/cluster%04d.svg", 
                  out.dir, k), toolTipMode = 2, title = sprintf("Bicluster %04d", 
                  k), xmlHeader = T)
                plotClust(k, w.motifs = T, seq.type = seq.type, 
                  ...)
                dev.off()
            })
        }
        cat("\n")
    }
    if ("pdf" %in% output) {
        require(igraph)
        cat("PDFS: ")
        lapply(ks, function(k) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/pdfs/cluster%04d.pdf", 
                out.dir, k))) 
                return(NULL)
            if (has.pdftk) 
                pdf(sprintf("%s/pdfs/cluster%04d.pdf", out.dir, 
                  k))
            else cairo_pdf(sprintf("%s/pdfs/cluster%04d.pdf", 
                out.dir, k))
            try(plotClust(k, w.motifs = T, seq.type = seq.type, 
                ...), silent = T)
            dev.off()
        })
        cat("\n")
    }
    if (gaggle && "html" %in% output) {
        require(hwriter)
        cat("HTMLS: ")
        lapply(ks, function(k, ...) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/htmls/cluster%04d.html", 
                out.dir, k))) 
                return()
            rows <- sort(get.rows(k))
            if (length(rows) <= 0) 
                return()
            short.names <- get.long.names(rows, short = T)
            short.names <- cbind(rows, short.names)
            rownames(short.names) <- colnames(short.names) <- NULL
            long.names <- get.long.names(rows, short = F)
            long.names <- cbind(rows, long.names)
            rownames(long.names) <- colnames(long.names) <- NULL
            refseq.names <- unique(unlist(get.synonyms(rows)))
            refseq.names <- grep("^NP_", refseq.names, val = T)
            upstream.seqs <- try(get.sequences(k, filter = F, 
                uniq = F), silent = T)
            if (class(upstream.seqs) == "try-error" || is.null(upstream.seqs)) {
                upstream.seqs <- rep("", length(rows))
                names(upstream.seqs) <- rows
            }
            upstream.seqs <- cbind(names(upstream.seqs), upstream.seqs)
            rownames(upstream.seqs) <- colnames(upstream.seqs) <- NULL
            htmltext <- paste(c("<html><head><title>Bicluster %K (%FILE)</title>", 
                "<style type=\"text/css\">", "  .hidden {", "     display: none;", 
                "   }", "  .gaggle-data {", "     color: green;", 
                "     font-size: xx-small;", "   }", "   p {", 
                "     color: red;", "     font-size: x-small;", 
                "   }", "</style>", "<script type=\"text/javascript\">", 
                "   function toggleVisible(id){", "      if (document.getElementById){", 
                "         obj = document.getElementById(id);", 
                "         if (obj) {", "            if (obj.style.display == 'none'){", 
                "               obj.style.display = 'block';", 
                "            } else {", "               obj.style.display = 'none';", 
                "            }", "         }", "      }", "   }", 
                "</script>", "</head>", "<table><tr><td>", "<iframe src=\"../svgs/cluster%K03d%K.svg\" width=\"600\" height=\"520\" frameborder=\"0\"></iframe>", 
                "</td><td>", "<p><a href=\"#bicluster%K03d%K\" onclick=\"toggleVisible('bicluster%K03d%K'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows and columns.</p>", 
                "<div id=\"bicluster%K03d%K\" style=\"display:none;\" class=\"gaggle-data bicluster\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                  length(rows), length(get.cols(k))), "   <div class=\"gaggle-cluster\">", 
                "      <ol class=\"gaggle-rowNames\">", paste("<li>", 
                  sort(rows), "</li>", sep = "", collapse = ""), 
                "      </ol>", "   <ol class=\"gaggle-columnNames\">", 
                paste("<li>", sort(get.cols(k)), "</li>", sep = "", 
                  collapse = ""), "      </ol>", "   </div>", 
                "</div>", "<p><a href=\"#bicluster%K03d%K_genes\" onclick=\"toggleVisible('bicluster%K03d%K_genes'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (genes).</p>", 
                "<div id=\"bicluster%K03d%K_genes\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K genes</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  length(rows)), "   <div class=\"gaggle-namelist\">", 
                "      <ol>", paste("<li>", sort(rows), "</li>", 
                  sep = "", collapse = ""), "      </ol>", "   </div>", 
                "</div>", "<p><a href=\"#bicluster%K03d%K_short_names\" onclick=\"toggleVisible('bicluster%K03d%K_short_names'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (short gene names).</p>", 
                "<div id=\"bicluster%K03d%K_short_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K short names</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  nrow(short.names)), "   <span class=\"gaggle-namelist-tag hidden\">short_name</span>", 
                hwrite(short.names, table.class = "toc", col.class = list(NA, 
                  "short_name"), border = 1, table.style = "font-family: monospace; font-size: xx-small; color: green; border-collapse: collapse"), 
                "   </div>", "<p><a href=\"#bicluster%K03d%K_long_names\" onclick=\"toggleVisible('bicluster%K03d%K_long_names'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (long gene names).</p>", 
                "<div id=\"bicluster%K03d%K_long_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K long names</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  nrow(long.names)), "   <span class=\"gaggle-namelist-tag hidden\">long_name</span>", 
                hwrite(long.names, table.class = "toc", col.class = list(NA, 
                  "long_name"), border = 1, table.style = "font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse"), 
                "   </div>", "<p><a href=\"#bicluster%K03d%K_refseq_names\" onclick=\"toggleVisible('bicluster%K03d%K_refseq_names'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K rows (NCBI RefSeq gene IDs).</p>", 
                "<div id=\"bicluster%K03d%K_refseq_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K NCBI RefSeq IDs</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  length(refseq.names)), "   <div class=\"gaggle-namelist\">", 
                "      <ol>", paste("<li>", sort(refseq.names), 
                  "</li>", sep = "", collapse = ""), "      </ol>", 
                "   </div>", "</div>", "<p><a href=\"#bicluster%K03d%K_upstream_seqs\" onclick=\"toggleVisible('bicluster%K03d%K_upstream_seqs'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K gene upstream sequences.</p>", 
                "<div id=\"bicluster%K03d%K_upstream_seqs\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K upstream sequences</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  nrow(upstream.seqs)), "   <span class=\"gaggle-namelist-tag hidden\">upstream</span>", 
                hwrite(upstream.seqs, table.class = "toc", col.class = list(NA, 
                  "upstream"), border = 1, table.style = "font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse"), 
                "   </div>", "<p><a href=\"#bicluster%K03d%K_arrays\" onclick=\"toggleVisible('bicluster%K03d%K_arrays'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K columns (arrays; conditions).</p>", 
                "<div id=\"bicluster%K03d%K_arrays\" style=\"display:none;\" class=\"gaggle-data arrays\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K arrays</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%d</span>", 
                  length(get.cols(k))), "   <div class=\"gaggle-namelist\">", 
                "      <ol>", paste("<li>", sort(get.cols(k)), 
                  "</li>", sep = "", collapse = ""), "      </ol>", 
                "   </div>", "</div>", "<p><a href=\"#bicluster%K03d%K_ratios\" onclick=\"toggleVisible('bicluster%K03d%K_ratios'); return false;\">[+]</a>", 
                "Show/hide bicluster #%K ratios.</p>", "<div id=\"bicluster%K03d%K_ratios\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                "   <span class=\"gaggle-name hidden\">bicluster %K ratios</span>", 
                "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                  length(rows), length(get.cols(k))), "   <div class=\"gaggle-matrix-tsv\">", 
                "        RATIOS", "   </div>", "</div>", if (!is.null(seq.type) && 
                  !is.null(meme.scores[[seq.type]][[k]]$meme.out) && 
                  !is.null(meme.scores[[seq.type]][[k]]$meme.out[[1]])) paste("<p><a href=\"#bicluster%K03d%K_pssm1\" onclick=\"toggleVisible('bicluster%K03d%K_pssm1'); return false;\">[+]</a>", 
                  "Show/hide bicluster #%K motif PSSM #1.</p>", 
                  "<div id=\"bicluster%K03d%K_pssm1\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                  "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #1</span>", 
                  "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                  sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                    nrow(meme.scores[[seq.type]][[k]]$meme.out[[1]]$pssm), 
                    ncol(meme.scores[[seq.type]][[k]]$meme.out[[1]]$pssm)), 
                  "   <div class=\"gaggle-matrix-tsv\">", "           MOTIF1", 
                  "   </div>", "</div>") else "", if (!is.null(seq.type) && 
                  length(meme.scores[[seq.type]][[k]]$meme.out) >= 
                    2 && !is.null(meme.scores[[seq.type]][[k]]$meme.out[[2]])) paste("<p><a href=\"#bicluster%K03d%K_pssm2\" onclick=\"toggleVisible('bicluster%K03d%K_pssm2'); return false;\">[+]</a>", 
                  "Show/hide bicluster #%K motif PSSM #2.</p>", 
                  "<div id=\"bicluster%K03d%K_pssm2\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                  "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #2</span>", 
                  "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                  sprintf("   <span class=\"gaggle-size hidden\">%dx%d</span>", 
                    nrow(meme.scores[[seq.type]][[k]]$meme.out[[2]]$pssm), 
                    ncol(meme.scores[[seq.type]][[k]]$meme.out[[2]]$pssm)), 
                  "   <div class=\"gaggle-matrix-tsv\">", "           MOTIF2", 
                  "   </div>", "</div>") else "", "</td></table>", 
                if ("pdf" %in% output) sprintf("<a href=\"../pdfs/cluster%04d.pdf\">View PDF version</a>", 
                  k) else "", "</html>"), collapse = "\n")
            rm(short.names, long.names, refseq.names, upstream.seqs)
            htmltext <- gsub("%K03d%K", sprintf("%04d", k), htmltext)
            htmltext <- gsub("%K", k, htmltext)
            htmltext <- gsub("%FILE", cmonkey.filename, htmltext)
            htmltext <- gsub("%SPECIES", gsub("_", " ", rsat.species), 
                htmltext)
            tmp <- as.data.frame(get.cluster.matrix(rows, get.cols(k)))
            tmp <- cbind(GENES = rownames(tmp), tmp)
            tf <- tempfile()
            write.table(tmp, file = tf, sep = "\t", quote = F, 
                row.names = F)
            rm(tmp)
            htmltext <- sub("RATIOS", paste(readLines(tf), collapse = "\n"), 
                htmltext)
            unlink(tf)
            if (!is.null(seq.type) && !is.null(meme.scores[[seq.type]][[k]]$meme.out)) {
                if (!is.null(meme.scores[[seq.type]][[k]]$meme.out[[1]])) {
                  tmp <- as.data.frame(meme.scores[[seq.type]][[k]]$meme.out[[1]]$pssm)
                  if (!is.null(tmp) && nrow(tmp) > 0) {
                    tmp <- cbind(1:nrow(tmp), tmp)
                    colnames(tmp) <- c("POSITION", "A", "C", 
                      "G", "T")
                    write.table(tmp, file = tf, sep = "\t", quote = F, 
                      row.names = F)
                    htmltext <- sub("MOTIF1", paste(readLines(tf), 
                      collapse = "\n"), htmltext)
                    unlink(tf)
                  }
                  rm(tmp)
                }
                if (length(meme.scores[[seq.type]][[k]]$meme.out) >= 
                  2 && !is.null(meme.scores[[seq.type]][[k]]$meme.out[[2]])) {
                  tmp <- as.data.frame(meme.scores[[seq.type]][[k]]$meme.out[[2]]$pssm)
                  if (!is.null(tmp) && nrow(tmp) > 0) {
                    tmp <- cbind(1:nrow(tmp), tmp)
                    colnames(tmp) <- c("POSITION", "A", "C", 
                      "G", "T")
                    write.table(tmp, file = tf, sep = "\t", quote = F, 
                      row.names = F)
                    htmltext <- sub("MOTIF2", paste(readLines(tf), 
                      collapse = "\n"), htmltext)
                    unlink(tf)
                  }
                  rm(tmp)
                }
            }
            rm(tf)
            cat(htmltext, file = sprintf("%s/htmls/cluster%04d.html", 
                out.dir, k), sep = "\n")
            rm(htmltext)
        })
        cat("\n")
    }
    if ("png" %in% output) {
        mc <- get.parallel(length(ks), para = 1)
        cat("PROFILES: ")
        lapply(ks, function(k, ...) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/htmls/cluster%04d_profile.png", 
                out.dir, k))) 
                return()
            try({
                c <- get.clust(k)
                png(sprintf("%s/htmls/cluster%04d_profile.png", 
                  out.dir, k), width = 128, height = 64, antialias = "subpixel")
                par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 0) * 
                  0.75)
                plotCluster(c, main = "", no.par = T, ...)
                dev.off()
            }, silent = T)
        })
        cat("\n")
        cat("NETWORKS: ")
        require(igraph)
        lapply(ks, function(k, ...) {
            if (k%%25 == 0) 
                cat(k)
            else cat(".")
            if (file.exists(sprintf("%s/htmls/cluster%04d_network.png", 
                out.dir, k))) 
                return()
            try({
                png(sprintf("%s/htmls/cluster%04d_network.png", 
                  out.dir, k), width = 64, height = 64, antialias = "subpixel")
                par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 0) * 
                  0.75)
                c <- get.clust(k)
                plotCluster.network(c, cex = 0.3, no.legend = T, 
                  ...)
                dev.off()
            }, silent = T)
        })
        cat("\n")
        if (!is.null(seq.type)) {
            cat("MOTIFS: ")
            lapply(ks, function(k, ...) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                e.vals <- lapply(meme.scores[[seq.type]][[k]]$meme.out, 
                  "[[", "e.value")
                pssms <- lapply(meme.scores[[seq.type]][[k]]$meme.out, 
                  "[[", "pssm")
                if (length(pssms) < 2) {
                  for (i in (length(pssms) + 1):2) {
                    pssms[[i]] <- matrix(0.25, nrow = 6, ncol = 4)
                    e.vals[[i]] <- Inf
                  }
                }
                for (pp in 1:length(pssms)) {
                  if (file.exists(sprintf("%s/htmls/cluster%04d_pssm%d.png", 
                    out.dir, k, pp))) 
                    return()
                  try({
                    png(sprintf("%s/htmls/cluster%04d_pssm%d.png", 
                      out.dir, k, pp), width = 128, height = 64, 
                      antialias = "subpixel")
                    if (is.matrix(pssms[[pp]])) 
                      try(viewPssm(pssms[[pp]], e.val = NA, mot.ind = pp, 
                        main.title = sprintf("e=%.3g", e.vals[[pp]]), 
                        cex.main = 0.7), silent = T)
                    dev.off()
                  }, silent = T)
                }
            })
            cat("\n")
            cat("MOTIF POSITIONS: ")
            lapply(ks, function(k, ...) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (file.exists(sprintf("%s/htmls/cluster%04d_mot_posns.png", 
                  out.dir, k))) 
                  return()
                try({
                  png(sprintf("%s/htmls/cluster%04d_mot_posns.png", 
                    out.dir, k), width = 128, height = 12 + 6 * 
                    length(get.rows(k)), antialias = "subpixel")
                  par(mar = rep(0.5, 4) + 0.1, mgp = c(3, 1, 
                    0) * 0.75)
                  c <- plotClust(k, dont.plot = T, ...)
                  plotClusterMotifPositions(c, cex = 0.4, no.key = T, 
                    ...)
                  dev.off()
                }, silent = F)
            })
            cat("\n")
        }
    }
    if ("main" %in% output) {
        mc <- get.parallel(length(ks))
        cat("WRITING MAIN HTML TABLE...")
        require(hwriter)
        dlf(paste(out.dir, "hwriter.css", sep = "/"), "http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css")
        dlf(paste(out.dir, "sorttable.js", sep = "/"), "http://www.kryogenix.org/code/browser/sorttable/sorttable.js")
        cat("...")
        cluster.summ <- cluster.summary(e.cutoff = NA, nrow.cutoff = 2)
        write.table(cluster.summ, file = paste(out.dir, "/cluster.summary.tsv", 
            sep = ""), quote = F, sep = "\t")
        cat("...")
        html <- openPage(paste(out.dir, "/index.html", sep = ""), 
            link.javascript = "sorttable.js", title = paste("cMonkey bicluster summary for run", 
                cmonkey.filename), link.css = "hwriter.css")
        hwrite(paste("<h2>cMonkey bicluster summary for run", 
            cmonkey.filename, "</h2>"), html)
        hwrite("<ul><li>Download a tab-delimited version of this table", 
            html, link = "cluster.summary.tsv", style = "font-size:75%")
        hwrite("<li>Download a list of each bicluster's gene members", 
            html, link = "cluster.members.genes.txt", style = "font-size:75%")
        hwrite("<li>Download a list of each bicluster's array/condition members", 
            html, link = "cluster.members.arrays.txt", style = "font-size:75%")
        cat("...")
        hwrite("<li>Plots of summary statistics of biclustering run", 
            html, link = "svgs/stats.svg", style = "font-size:75%")
        hwrite("<li>Saved cMonkey R session file", html, link = "cm_session.RData", 
            style = "font-size:75%")
        hwrite("<li>Summary of cMonkey input parameters</ul>", 
            html, link = "cm.params.txt", style = "font-size:75%")
        hwrite("<br><center><b>Bicluster summary</b></center><br>", 
            html)
        hwrite("<br><center><b>Sort the table by a given column by clicking on the column's header.<br>Click on bicluster link in first column for more info.</b></center><br>", 
            html, style = "font-size:60%")
        cat("...")
        himg0 <- hwriteImage(sprintf("htmls/cluster%04d_profile.png", 
            as.integer(rownames(cluster.summ))), table = F)
        himg0 <- hwrite(paste(himg0, sprintf("Residual = %.3f", 
            cluster.summ$resid), sep = "<br>"), center = TRUE, 
            table = F)
        himg0a <- hwriteImage(sprintf("htmls/cluster%04d_network.png", 
            as.integer(rownames(cluster.summ))), table = F)
        if (!is.null(seq.type)) {
            e.val.1 <- lapply(meme.scores[[seq.type]][as.integer(rownames(cluster.summ))], 
                function(i) i$meme.out[[1]]$e.value)
            for (i in 1:length(e.val.1)) if (is.null(e.val.1[[i]])) 
                e.val.1[[i]] <- NA
            himg1 <- hwriteImage(sprintf("htmls/cluster%04d_pssm1.png", 
                as.integer(rownames(cluster.summ))), table = F, 
                title = sprintf("E-val = %.3g", unlist(e.val.1)))
            himg1 <- hwrite(paste(himg1, as.character(cluster.summ$consensus1), 
                sep = "<br>"), center = TRUE, table = F)
            if (!is.null(seq.type)) 
                e.val.2 <- lapply(meme.scores[[seq.type]][as.integer(rownames(cluster.summ))], 
                  function(i) i$meme.out[[2]]$e.value)
            else e.val.2 <- as.list(rep(NA, k.clust))
            for (i in 1:length(e.val.2)) if (is.null(e.val.2[[i]])) 
                e.val.2[[i]] <- NA
            himg2 <- hwriteImage(sprintf("htmls/cluster%04d_pssm2.png", 
                as.integer(rownames(cluster.summ))), table = F, 
                title = sprintf("E-val = %.3g", unlist(e.val.2)))
            himg2 <- hwrite(paste(himg2, as.character(cluster.summ$consensus2), 
                sep = "<br>"), center = TRUE, table = F)
            himg2a <- hwriteImage(sprintf("htmls/cluster%04d_mot_posns.png", 
                as.integer(rownames(cluster.summ))), table = F)
            e.val.1[is.na(e.val.1)] <- 9e+09
            e.val.2[is.na(e.val.2)] <- 9e+09
        }
        else {
            e.val.1 <- e.val.2 <- as.list(rep(NA, k.clust))
            himg1 <- himg2 <- himg2a <- NULL
        }
        cluster.summ$score <- sprintf("%.3f", cluster.summ$score)
        rn <- rownames(cluster.summ)
        cat("...")
        cluster.summ.orig <- cluster.summ
        cluster.summ <- cbind(bicluster = cluster.summ$k, n.genes = cluster.summ$nrow, 
            n.arrays = sapply(as.integer(rownames(cluster.summ)), 
                function(i) length(get.cols(i))), score = cluster.summ$score, 
            residual = sprintf("%.3f", cluster.summ$resid))
        if ("score.norm" %in% colnames(cluster.summ.orig)) 
            cluster.summ <- cbind(cluster.summ, score.norm = sprintf("%.3f", 
                cluster.summ.orig$score.norm))
        rownames(cluster.summ) <- rn
        rows <- list()
        for (k in as.integer(rn)) rows[[k]] <- sort(get.rows(k))
        himg3 <- hwrite(sapply(as.integer(rn), function(k) paste(rows[[k]], 
            collapse = " ")), table = F)
        cat("...\n")
        if (!no.genome.info) {
            himg4 <- hwrite(unlist(mc$apply(as.integer(rn), function(k) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (length(rows[[k]]) <= 0) 
                  return()
                tmp <- get.long.names(rows[[k]], short = T)
                tmp <- unique(tmp[!tmp %in% rows[[k]] & tmp != 
                  ""])
                paste(tmp, collapse = " ")
            })), table = F)
            cat("\n")
            himg5 <- hwrite(unlist(mc$apply(as.integer(rn), function(k) {
                if (k%%25 == 0) 
                  cat(k)
                else cat(".")
                if (length(rows[[k]]) <= 0) 
                  return()
                tmp <- get.long.names(rows[[k]], short = F)
                tmp <- unique(tmp[!tmp %in% rows[[k]] & tmp != 
                  ""])
                paste(tmp, collapse = " | ")
            })), table = F)
            cat("\n")
        }
        else {
            himg4 <- himg5 <- NULL
        }
        nas <- rep(NA, nrow(cluster.summ))
        hwrite(cbind(cluster.summ[, 1:min(ncol(cluster.summ), 
            6)], profile = himg0, network = himg0a, motif1 = himg1, 
            motif2 = himg2, motif.posns = himg2a, probe.names = himg3, 
            short.names = himg4, long.names = himg5), html, row.names = F, 
            table.style = "text-align:center;font-size:70%;font-family:Arial", 
            table.class = "sortable", row.style = list("font-weight:bold;text-align:center;font-size:70"), 
            col.style = list(probe.names = "font-size:70%", orf.names = "font-size:50%", 
                short.names = "font-size:50%", long.names = "font-size:50%", 
                motif1 = "font-size:50%", motif2 = "font-size:50%"), 
            col.sorttable_customkey = list(residual = sprintf("%.3f", 
                cluster.summ.orig$residual), score.norm = if ("score.norm" %in% 
                colnames(cluster.summ.orig)) sprintf("%.3f", 
                cluster.summ.orig$score.norm) else NULL, profile = sprintf("%.3f", 
                cluster.summ.orig$resid), motif1 = sprintf("%.30f", 
                unlist(e.val.1)), e.val1 = sprintf("%.30f", unlist(e.val.1)), 
                motif2 = sprintf("%.30f", unlist(e.val.2)), e.val2 = sprintf("%.30f", 
                  unlist(e.val.2))), col.class = list(network = c("sorttable_nosort", 
                nas), motif.posns = c("sorttable_nosort", nas)), 
            col.link = list(sprintf("htmls/cluster%04d.html", 
                as.integer(rownames(cluster.summ)))))
        closePage(html, splash = F)
        for (i in sapply(1:k.clust, function(k) c(k, sort(get.rows(k))))) cat(i, 
            "\n", file = paste(out.dir, "/cluster.members.genes.txt", 
                sep = ""), append = T)
        for (i in sapply(1:k.clust, function(k) c(k, sort(get.cols(k))))) cat(i, 
            "\n", file = paste(out.dir, "/cluster.members.arrays.txt", 
                sep = ""), append = T)
        tmp <- capture.output(for (name in ls(cmonkey.params)) {
            cat(name, "= ")
            str(get(name, envir = cmonkey.params), no.list = T)
        })
        cat(tmp, file = paste(out.dir, "/cm.params.txt", sep = ""), 
            sep = "\n", collapse = "\n")
    }
    if (gzip) {
        rpl <- function(find, replace, file, ...) {
            f <- readLines(file)
            f <- gsub(find, replace, f, ...)
            writeLines(f, con = file)
        }
        system(sprintf("gzip -v %s/svgs/*.svg", out.dir))
        for (f in list.files(paste(out.dir, "/svgs", sep = ""), 
            full = T)) if (grepl(".svg.gz", f, fixed = T)) 
            system(sprintf("mv -v %s %s", f, sub(".svg.gz", ".svgz", 
                f, fixed = T)))
        lapply(c(list.files(sprintf("%s/htmls", out.dir), pattern = glob2rx("*.html"), 
            full = T), list.files(out.dir, pattern = glob2rx("*.html"), 
            full = T)), function(f) {
            cat(f, "\n")
            rpl(".svg\"", ".svgz\"", f, fixed = T)
        })
        if (has.pdftk) 
            lapply(list.files(paste(out.dir, "/pdfs", sep = ""), 
                full = T), function(f) if (grepl(".pdf", f, fixed = T)) {
                system(sprintf("pdftk %s output %s.tmp compress", 
                  f, f))
                system(sprintf("/bin/mv -fv %s.tmp %s", f, f))
            })
    }
    if ("rdata" %in% output) 
        save.cmonkey.env(file = paste(out.dir, "/cm_session.RData", 
            sep = ""))
    out.dir
}
