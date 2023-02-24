bam.filter.report <- function(chunk.l, project.name, report.only = F, filter.col = NULL, accept.region = NULL,
                              remove.mate = T, out.dir, threads = 6,
                              header.file, tags, shift) {
  warning("this function assume that one read wouldn't align twice within each interval")
  if (remove.mate == T) {
    warning("using the primitive way to remove mate, which is to remove all alignments with the same qname")
  }
  if (is.data.frame(chunk.l)) {
    stop("must be a chunk list")
  }
  
  if (is.null(names(chunk.l))) {
    stop("chunk list must be named!!")
  }
  
  dir.create(out.dir, showWarnings = F, recursive = T)
  
  res <- mclapply(names(chunk.l), function(chunk.name) {
    chunk <- chunk.l[[chunk.name]]
    if (report.only == F) {
      chunk <- utilsFanc::add.column.fanc(df1 = chunk, df2 = data.frame(pass = rep(NA, nrow(chunk))), pos = 1)
      fail.qnames <- chunk$qname[chunk[, filter.col] < accept.region[1] | chunk[,filter.col] > accept.region[2]]
      chunk$pass[chunk$qname %in% fail.qnames] <- F
      chunk$pass[! chunk$qname %in% fail.qnames] <- T
    }
    
    stat <- chunk %>% group_by(pass) %>% 
      summarise(n.reads = n(), frac.secondary = sum(isSecondaryAlignment)/n(),
                frac.mapq0 = sum(min.mapq == 0)/n(), frac.mapq1 = sum(min.mapq == 1)/n(), 
                frac.mapqh = sum(min.mapq > 1)/n(), qwidth = mean(qwidth), 
                isize = mean(isize), m.score = mean(m.score)) %>% 
      as.data.frame() %>% 
      mutate_if(is.numeric, round, digits = 3)
    stat <- utilsFanc::add.column.fanc(df1 = stat, df2 = data.frame(name = rep(chunk.name, nrow(stat))),
                                       pos = 1)  
    # write.table(stat, paste0(out.dir, "/", project.name, "_", chunk.name, "_stat.tsv"), quote = F, col.names = T, 
    #             row.names = F, sep = "\t")
    res <- list(pass = chunk[chunk$pass == T, ], fail = chunk[chunk$pass == F, ], stat = stat,
                filter.col = filter.col, accept.region = accept.region)
    return(res)
  }, mc.cores = threads, mc.cleanup = T)
  names(res) <- names(chunk.l)
  
  saveRDS(res, paste0(out.dir, "/", project.name, ".Rds"))
  
  stat <- lapply(res, function(x) return(x$stat)) %>% Reduce(rbind, .)
  write.table(stat, paste0(out.dir, "/", project.name, "_stat.tsv"), quote = F, col.names = T, 
              row.names = F, sep = "\t")

  to.bam <- list(pass = lapply(res, function(x) return(x$pass)),
                 fail = lapply(res, function(x) return(x$fail)))
  
  bams <- sapply(names(to.bam), function(x) {
    bam <- write.bam(bam.df = to.bam[[x]], out.bam = paste0(out.dir, "/", project.name, "_", x, ".bam"),
              header.file = header.file, fields = ALL.FIELDS,
              tags = tags, append = F, sam2bam = T, sort = T, index = T,
              shift = shift)
    return(bam)
  })
  cat(utilsFanc::bash2ftp(bams))
  return(res)
}

nkc.composite.filter <- function(chunk.l, m.score.cutoff, out.dir, 
                                 header.file, shift, tags, threads = 1) {
  # first filter by m.score
  chunk.l <- bam.filter.report(chunk.l = chunk.l, project.name = "m.score_filter", 
                               report.only = F, filter.col = "m.score", 
                               accept.region = c(m.score.cutoff, 1), 
                               remove.mate = T, out.dir = out.dir,
                               tags = tags, shift = shift, header.file = header.file,
                               threads = threads) %>% 
    lapply(function(x) return(x$pass))
  
  chunk.l <- mclapply(chunk.l, function(chunk) {
    chunk$pass <- NA
    chunk$pass[chunk$isSecondaryAlignment == F] <- T
    chunk$pass[chunk$isSecondaryAlignment == T & chunk$min.mapq < 2] <- T
    chunk$pass[is.na(chunk$pass)] <- F
    chunk$pass[chunk$tag.YT != "CP"] <- F
    return(chunk)
  }, mc.cores = threads, mc.cleanup = T)
  res <-  bam.filter.report(chunk.l = chunk.l, project.name = "second_filter", 
                            report.only = T, # filter.col = "m.score", 
                            # accept.region = c(m.score.cutoff, 1), 
                            remove.mate = T, out.dir = out.dir,
                            tags = tags, shift = shift, header.file = header.file,
                            threads = threads)
  return(res)
}