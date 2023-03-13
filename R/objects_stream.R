bam.subset.wrapper <- function(bam.vec, out.bam.vec = NULL, out.dir, pct,
                               threads.bam = 4, threads.sub = 4,
                               samtools = SAMTOOLS, seed= 42, run = T) {
  if (is.null(out.bam.vec)) {
    out.bam.vec <- paste0(out.dir, "/", basename(bam.vec) %>% sub(".bam$", paste0(".subset",pct,".bam"), .))
  }
  pct <- as.character(pct) %>% sub("0.", "", .)
  out.bams <- mclapply(seq_along(bam.vec), function(i) {

    bam <- bam.vec[i]
    out.bam <- out.bam.vec[i]
    dir.create(dirname(out.bam), showWarnings = F, recursive = T)
    cmd <- paste0(samtools, " view -h -b -o ", out.bam, " -@ ", threads.sub,
                  " -s ", seed, ".", pct, " ", bam)
    utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = run)
    if (!file.exists(out.bam))
      stop(paste0(out.bam, " failed to generate"))
    return(out.bam)
  }, mc.cleanup = T, mc.cores = threads.bam)
  return(out.bams)
}


bam.chunk.to.df <- function(chunk) {
  null.names <- names(chunk[[1]]$tag[sapply(chunk[[1]]$tag, is.null)])
  for (i in null.names) {
    chunk[[1]]$tag[[i]] <- rep(NA, length(chunk[[1]]$qname))
  }
  chunk <- as.data.frame(chunk)
  return(chunk)
}

bam.chunk.to.df.2 <- function(chunks, rbind = T,
                              add.bits = NULL,
                              shift = 0, shift.cols = c("pos", "mpos"),
                              add.shift = 0, add.min.mapq = F, add.m.score = F) {
  res <- lapply(chunks, function(chunk) {
    chunk <- scFanc::factor2character.fanc(chunk)
    null.names <- names(chunk$tag[sapply(chunk$tag, is.null)])
    for (i in null.names) {
      chunk$tag[[i]] <- rep(NA, length(chunk$qname))
    }
    chunk <- as.data.frame(chunk)
    if (!is.null(add.bits))
      chunk <- bam.df.add.flag(bam.df = chunk, flag.bits = add.bits, pos = 1)
    if (add.min.mapq == T) {
      chunk <- chunk %>% arrange(qname) %>% group_by(qname) %>% mutate(min.mapq = min(mapq)) %>%
        ungroup() %>% as.data.frame()
      chunk <- utilsFanc::df.rearrange.cols(df = chunk, cols = "min.mapq", after = "mapq")
    }
    if (shift != 0)
      chunk <- bam.shift(bam.df = chunk, shift = shift, shift.cols = shift.cols)
    if (add.shift != 0) {
      chunk <- utilsFanc::add.column.fanc(df1 = chunk, df2 = data.frame(pos.shift = chunk$pos + add.shift,
                                                                            end = chunk$pos + add.shift + 1000),
                                            after = "pos")
    }
    if (add.m.score == T) {
      chunk <- utilsFanc::add.column.fanc(df1 = chunk, df2 = data.frame(qwidth = chunk$seq %>% nchar(),
                                                                            m.score = chunk$tag.AS/(chunk$seq %>% nchar())),
                                            after = "min.mapq")
    }
    return(chunk)
  })
  if (rbind == T) {
    res <- res %>% Reduce(rbind,.) %>% unique()
  }
  return(res)
}

bam.df.add.flag <- function(bam.df, flag.bits, ...) {
  df2 <- Rsamtools::bamFlagAsBitMatrix(flag = bam.df$flag %>% as.integer())[, flag.bits] %>%
    as.data.frame() %>% `colnames<-`(flag.bits)
  bam.df <- utilsFanc::add.column.fanc(df1 = bam.df, df2 = df2, ...)
  return(bam.df)
}

bam.df.subset.by.region <- function(bam.df, gr, return.list = F, gr.name.col = NULL,
                                    rescue.qname = T) {
  tmp <- bam.df
  tmp$pos.end <- tmp$pos + GenomicAlignments::cigarWidthAlongReferenceSpace(tmp$cigar)
  tmp <- tmp %>% dplyr::select(rname, pos, pos.end, qname)
  colnames(tmp) <- c("chr", "start", "end", "qname")
  tmp <- makeGRangesFromDataFrame(tmp, keep.extra.columns = T)
  o <- findOverlaps(tmp, gr)
  if (is.null(gr.name.col)) {
    gr$locus <- utilsFanc::gr.get.loci(gr)
    gr.name.col <- "locus"
  }
  bam.list <- queryHits(o) %>% split(f = mcols(gr)[, gr.name.col][subjectHits(o)]) %>%
    lapply(function(x) {
      if (rescue.qname == T) {
        qnames.in <- tmp[x]$qname
        bam.df <- bam.df %>% filter(qname %in% qnames.in)
      } else {
        bam.df <- bam.df[unique(x),]
      }
      return(bam.df)
    })


  if (length(bam.list) == 1)
    bam.list <- bam.list[[1]]
  return(bam.list)

  # if (return.list == F) {
  #   qnames.in <- tmp[queryHits(o) %>% unique()]$qname
  #   bam.df <- bam.df %>% filter(qname %in% qnames.in)
  #   return(bam.df)
  # } else {

  # }
}

write.bam <- function(bam.df, out.bam, header.file, fields = ALL.FIELDS, tags,
                      shift = 0,
                      append = F, sam2bam = T,
                      sort = T, index = T,
                      threads = 4, samtools = SAMTOOLS) {
  if (shift != 0) {
    warning("when shifting reads, make sure to use the correct header that will correspond to the shifted coordinates!!")
  }
  if (is.list(bam.df) && !is.data.frame(bam.df))
    bam.df <- Reduce(rbind, bam.df)
  bam.df <- bam.df[, c(fields, paste0("tag.", tags))] %>%
    tag.formatting() %>% rnext.formatting() %>% bam.shift(shift)
  out.sam <- sub(".bam$", ".sam", out.bam)
  if (append == F) {
    cmd <- paste0("rsync ", header.file, " ", out.sam)
    system(cmd)
  }
  write.table(bam.df, file = out.sam, append = T, quote = F, sep = "\t", row.names = F, col.names = F)
  if (sam2bam == T) {
    if (sort == T) {
      cmd <- paste0(samtools," sort -m 2G -@ ", threads, " ", out.sam, " > ", out.bam)
    } else {
      cmd <- paste0(samtools," view -b -h -@ ", threads, " ", out.sam, " > ", out.bam)
    }
    system(cmd)
    if (index == T) {
      system(paste0(samtools, " index ", out.bam))
    }
  }
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " failed to generate"))
  else
    return(out.bam)
}


split.bam <- function(bam,  split.dir, split.size, override = F, run = T,
                      samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  prefix <- basename(bam) %>% paste0("_split_")
  system(paste0("mkdir -p ", split.dir))
  is.exist <- Sys.glob(paste0(split.dir, "/", prefix, "*"))
  if (length(is.exist) > 0) {
    if (override != T)
      stop(paste0("It seems that files with provided split.dir and prefix already exist\n",
                  "To avoid errors with globing later, the processed is stopped.\n",
                  "These files are: \n", paste0(is.exist, collapse = "\n")))
    else {
      cat(paste0("It seems that files with provided split.dir and prefix already exist\n",
                  "To avoid errors with globing later, these files will be removed pre flight.\n",
                  "These files are: ", paste0(is.exist, collapse = "\n"), "\n"))
      rm.cmd <- paste0("rm -rf ", paste0(is.exist, collapse = " "))
      print(Sys.time())
      print(rm.cmd)
      if (run == T)
        system(rm.cmd)
    }
  }

  cmd.1 <- paste0("cd ", split.dir,
                  " && ",samtools," view -H ", bam," > header")
  print(Sys.time())
  print(cmd.1)
  if (run == T)
    system(cmd.1)
  cmd.2 <- paste0("cd ", split.dir,
                  " && ",samtools," view ", bam, " | split - ",prefix,
                  " -l ",split.size ," --filter='cat header - | ",samtools," view -b -h - > $FILE.bam'")
  print(Sys.time())
  print(cmd.2)
  if (run == T)
    system(cmd.2)
  splits <- Sys.glob(paste0(split.dir, "/", prefix, "*"))
  return(splits)
}

sort.bam <- function(bam, outfile, is.name.sort, thread, tag = NULL, mem.per.thread,run = T,
                     samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  cmd <- paste0(samtools, " sort ", " -@ ", thread, " -m ", mem.per.thread)
  if (is.name.sort == T)
    cmd <- paste0(cmd, " -n ")
  if (!is.null(tag))
    cmd <- paste0(cmd, " -t ", tag)
  cmd <- paste0(cmd, " ", bam,  " > ", outfile)
  print(Sys.time())
  print(cmd)
  if (run == T)
    system(cmd)
  if (file.exists(outfile))
    return(outfile)
  else
    stop(paste0(Sys.time(), "\n" ,bam, " was not successfully sorted. ",
                outfile, " was not generated"))

}

bam.shift <- function(bam.df, shift = 0, shift.cols = c("pos", "mpos")) {
  if (shift == 0)
    return(bam.df)
  for (field in shift.cols) {
    bam.df[, field] <- bam.df[, field] + shift
  }
  return(bam.df)
}


stream.bam.core <- function(bam.file, fields, tags= character(length = 0),
                            header.file = NULL,
                            samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools",
                            chunk.size,
                            chunk.call.back,
                            chunk.call.back.param.list) {
  # chunk.call.back must be able to receive hinge.pre as a prameter, process the chunk and merge it with
  ##hinge.pre, and return a hinge.pre for the next cycle
  # it must also accept a binary "eof" parameter, which is a flag indicating that the input file has ended and
  ##enable returning a value that will be the final output of the stream.
  # it must also accept a binary "bof" parameter, which is a flag indicating that this is the first chunk.
  ##The first chunk might be processed differently, such as creating/removing files
  if (is.null(header.file)) {
    header.file <- bam.get.header(bam = bam.file, out.header = tempfile(), samtools = samtools)
  }
  param <- ScanBamParam(what = fields, tag = tags)
  bf <- open(BamFile(bam.file, yieldSize = chunk.size))
  hinge.pre <- NULL
  bof <- T
  while (1) {
    chunk <- scanBam(bf, param = param)

    if (length(chunk[[1]]$qname) < 1) {
      final.out <- do.call(chunk.call.back,
                           c(chunk.call.back.param.list,
                             list(chunk = chunk, eof = T, bof=bof,
                                  hinge.pre = hinge.pre, header.file = header.file)))
      return(final.out)
    }

    null.names <- names(chunk[[1]]$tag[sapply(chunk[[1]]$tag, is.null)])
    for (i in null.names) {
      chunk[[1]]$tag[[i]] <- rep(NA, length(chunk[[1]]$qname))
    }
    chunk <- as.data.frame(chunk)



    hinge.pre <- do.call(chunk.call.back,
                         c(chunk.call.back.param.list,
                           list(chunk = chunk, eof = F, bof = bof,
                                hinge.pre = hinge.pre, header.file = header.file)))
    bof <- F
  }
  close(bf)
  return(final.out)
}



t.f <- function(chunk, eof, bof, hinge.pre, outfile) {
  # just write unique reference name. remove hinge reference name
  append.mode <- T
  if (bof == T)
    append.mode <- F
  rname <- chunk %>% as.data.frame() %>% pull(rname) %>% sort() %>%
  if (! is.null(hinge.pre)) {
    if (hinge.pre == rname[1])
      rname[-1]
  }


}

bam2bed.chunk <- function(chunk, eof, bof, hinge.pre, outfile) {
  # note: this function will group by BC and UB. it will take in hinge.pre as the first line.
  # however, it will not write down the last UB/CB ground and instead pass it down as return value (hinge.pre for the next cycle)
  # above 2 lines of comments obsolete.
  append.mode <- T
  if (bof == T)
    append.mode <- F
  # first perform some basic filtering
  df <- chunk %>% filter(!is.na(tag.UB), !is.na(tag.CB), is.na(tag.GX), is.na(tag.GN),
                         mapq == 255) %>%
    mutate(right = GenomicAlignments::cigarWidthAlongReferenceSpace(cigar = cigar),
           cb.ub = paste0(tag.CB, "|", tag.UB)) %>%
    select(rname, pos, right, cb.ub, flag, strand, qname)
  write.table(df, outfile, append = append.mode, quote = F, col.names = F,
              row.names = F, sep = "\t")
  return(NULL)

}

bam.filter.chunk.bk <- function(chunk, eof, bof, hinge.pre, outfile, header.file, thread.sub=4) {
  out.tempt <- paste0(outfile, ".sam")
  if (eof == T) {
    cmd <- paste0("samtools view -b -h -@ ", thread.sub, " > ", outfile)
    print(cmd)
    system(cmd)
    return(NULL)
  }
  append.mode <- T
  if (bof == T) {
    system(paste0("cp ", header.file, " ", out.tempt))
  }
  # first perform some basic filtering
  df <- chunk %>% filter(!is.na(tag.UB), !is.na(tag.CB), is.na(tag.GX), is.na(tag.GN),
                         #is.na(tag.AN), is.na(tag.TX),
                         mapq == 255)
  if (nrow(df) < 1)
    return(NULL)
  else {
    df <- df %>% mutate(tag.GX = NULL,
                        tag.AN = NULL,
                        tag.TX = NULL,
                        tag.GN = NULL) %>%
      tag.formatting() %>%
      mutate(rnext = "*", pnext = 0, tlen = 0, seq = seq, qual = qual) %>%
      select(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual,
             tag.CB, tag.UB, tag.BC, tag.QT, tag.CR, tag.UR, tag.CY, tag.UY)
    write.table(df, out.tempt, append = T, quote = F, col.names = F,
                row.names = F, sep = "\t")
  }

  return(NULL)
}

tag.formatting <- function(df) {
  cnames <- colnames(df)
  df$tags <- "miao"
  for (tag in cnames[grepl("^tag.", cnames)]) {
    na.idx <- is.na(df[, tag])
    type <- if_else(is.character(df[, tag]), "Z", "i")
    df[,tag] <- paste0(sub("^tag.", "",tag), ":", type, ":", df[,tag])
    df[, tag][na.idx] <- "miao"
    df$tags <- paste0(df$tags, "\t", df[, tag])
    df[, tag] <- NULL
  }
  # browser()
  df$tags <- df$tags %>% gsub("miao\t*", "", .) %>% sub("\t$", "",.)
  return(df)
}

rnext.formatting <- function(df) {
  df <- scFanc::factor2character.fanc(df)
  df$mrnm[df$mrnm == df$rname] <- "="
  return(df)
}

bam.filter.master <- function(bam, outfile, chunk.size = 100000) {
  system(paste0("mkdir -p ", dirname(outfile)))
  stream.bam.core(bam.file = bam,
                  fields =  c("qname", "flag", "rname", "pos","mapq", "cigar", "seq", "qual"),
                  tags = c("CB", "UB","GX", "GN", "AN", "TX") %>% c("CR", "CY", "UR", "UY", "BC", "QT", "CB") %>% unique(),
                  chunk.size = 100000,
                  chunk.call.back = bam.filter.chunk,
                  chunk.call.back.param.list = list(outfile = outfile))
}

bam.tag.detect <- function(bam, tempt.file = NULL, tag.txt = NULL,
                           samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  if (is.null(tempt.file))
    tempt.file <- tempfile()
  cmd <- paste0(samtools, " view ", bam, " | ",
                "cut -f12- | sed 's/\\([^\\t]*\\):[^\\t]*:[^\\t]*/\\1/g' > ", tempt.file)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)

  tags <- readLines(tempt.file)
  tags <- strsplit(tags, "\t") %>% unlist() %>% unique()
  if (!is.null(tag.txt)) {
    write(tags, tag.txt, sep = "\n")
  }
  return(tags)
}

bam.get.header <- function(bam, out.header,
                           samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  cmd <- paste0(samtools," view -H ", bam," > ", out.header)
  utilsFanc::cmd.exec.fanc(cmd, intern = F, run = T)
  return(out.header)
}

bam.2.fastq.wrapper <- function(bam.file.vec, out.dir, fastq.root.vec = NULL,
                                fastq.ext.vec, threads = 4, bedtools = BEDTOOLS,
                                run = T) {
  stop("this function wouldn't really work. You need to add name sort functionality")
  if (is.null(fastq.root.vec)) {
    fastq.root.vec <- paste0(out.dir, "/", sub(".bam$", "", basename(bam.file.vec)))
  }
  mclapply(seq_along(bam.file.vec), function(i) {
    bam <- bam.file.vec[i]
    fastq.root <- fastq.root.vec[i]
    dir.create(dirname(fastq.root)[1], showWarnings = F, recursive = T)
    read1 <- paste0(fastq.root, fastq.ext.vec[1])
    cmd <- paste0(bedtools, " bamtofastq -i ", bam, " -fq ", read1)
    if (!is.na(fastq.ext.vec[2])) {
      cmd <- paste0(cmd, " -fq2 ", fastq.root, fastq.ext.vec[2])
    }
    utilsFanc::cmd.exec.fanc(cmd = cmd, run = run, intern = F)
    if (!file.exists(read1)) {
      stop(paste0(read1, " failed to generate"))
    }
    return()
  }, mc.cleanup = T, mc.cores = threads)
}

bam.2.fastq <- function(bam.df, out.fastq.root, pe,
                        fastq.ext.names = c("_1.fastq.gz", "_2.fastq.gz")) {
  stop("this function doesn't work")
  if (is.list(bam.df) && !is.data.frame(bam.df))
    bam.df <- Reduce(rbind, bam.df)
  bam.df <- bam.df %>% select(qname, seq, qual, flag) %>%
    mutate(plus = "+", qname = paste0("@", qname))
  if (pe == T) {
    bam.df <- bamFanc::bam.df.add.flag(bam.df = bam.df, flag.bits = "isFirstMateRead")
    read1.df <- bam.df$isFirstMateRead == 1
    read2.df <- bam.df$isFirstMateRead == 0
    read1.qnames <- read1.df$qname
    read2.qnames <- read2.df$qname
    qnames <- intersect(read1.qnames, read2.qnames)
    read1.df <- read1.df %>% filter(qname %in% qnames) %>% arrange(qname)
    read2.df <- read2.df[qnames, ]
    df.l <- list(read1.df, read2.df)
    names(df.l) <- fastq.ext.names
    out.files <- mclapply(fastq.ext.names, function(ext.name) {
      df <- df.l[[ext.name]] %>% select(qname, seq, plus, qual)
      outfile <- paste0(out.fastq.root, ext.name)
      write.table(df, file = outfile, row.names = F, col.names = F, sep = "\n", quote = F)
      system(paste0("gzip ", outfile))
      return(outfile)
    }, mc.cores = 2, mc.cleanup = T) %>% unlist()

    if (any(!file.exists(out.files))) {
      stop(paste0(paste0(out.files, collapse = " or "), " failed to generate"))
    }
    return(out.files)
  }
}

bam.get.cut.site.ez <- function(chunk.l, is.open.bed = F,
                                overall.shift = 0, left.shift = -4, right.shift = 5,
                                gr.region.name.map = NULL, name.col = NULL) {
  # this function can also process the open.bed file generated by the target pipeline.
  switch <- F
  if (is.data.frame(chunk.l) | is(chunk.l, "GRanges")) {
    chunk.l <- list(chunk1 = chunk.l)
    switch <- T
  }
  chunk.names <- names(chunk.l)
  chunk.l <- lapply(chunk.names, function(region) {
    chunk <- chunk.l[[region]]
    if (is.open.bed == T) {
      if (length(chunk) == 0)
        chunk$name <- character(0)
      gr <- (shift(chunk, overall.shift) - 74.5) %>% plyranges::mutate(qname = name) %>% plyranges::select(qname)
      df <- as.data.frame(gr %>% `names<-`(NULL)) %>% dplyr::rename(chr = seqnames)
    } else {
      df <- chunk %>% dplyr::select(qname, isFirstMateRead, rname, pos, mpos, cigar)
      bLeft <- (df$pos < df$mpos) | (df$isFirstMateRead == 1 & df$pos == df$mpos)
      df$pos.use <- NA
      df$shift <- NA
      df[bLeft, ] <- df[bLeft, ] %>%
        mutate(pos.use = pos, shift = left.shift)
      df[!bLeft, ] <- df[!bLeft, ] %>%
        mutate(pos.use = pos + GenomicAlignments::cigarWidthAlongReferenceSpace(cigar), shift = right.shift)

      df <- df %>% mutate(start = pos.use + overall.shift + shift) %>% mutate(end = start) %>%
        dplyr::select(rname, start, end, qname) %>% dplyr::rename(chr = rname)
      gr <- df %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    }



    if (!is.null(gr.region.name.map)) {
      # browser()
      region.gr <- gr.region.name.map %>% plyranges::filter(!!as.name(name.col) == region)
      if (length(region.gr) != 1) {
        stop("error at gr.region.name.map")
      }
      o <- findOverlaps(gr, region.gr, ignore.strand = T)
      gr <- gr[queryHits(o), ]
      df <- df[queryHits(o), ]
    }
    df.pe <- df %>% group_by(qname) %>% summarise(start.left= min(start),
                                                  start.right = max(start), n.mates = n()) %>%
      ungroup() %>% as.data.frame()
    res <- list(gr = gr, df = df, df.pe = df.pe)
    return(res)
  })
  names(chunk.l) <- chunk.names
  if (switch == T)
    chunk.l <- chunk.l[[1]]
  return(chunk.l)
}

cut.site.pipe <- function(bam.files, sample.names = NULL,
                          ctrl.bed.files = NULL,
                          region.gr, name.col,  shift = 0, padding = 5,
                          report.dir,
                          tag, header.file, m.score.cutoff,
                          threads.bam = 4, threads.region = 4) {
  print("this function uses shift differently. " %>%
          paste0("It caters the fact that region.gr should be in genomic coordinates",
                 ", and the bam is on a mini genome"))
  utilsFanc::write.zip.fanc(df = region.gr, out.file = paste0(report.dir, "/region.bed"), bed.shift = T)
  if (is.null(sample.names))
    sample.names <- names(bam.files)
  if (is.null(sample.names))
    stop("samples must be named")
  names(bam.files) <- sample.names
  cutsite <- mclapply(sample.names, function(sample) {
    bam <- bam.files[sample]
    region.gr.shift <- shift(region.gr, -shift)
    # browser()
    seqlevels(region.gr) <- seqnames(region.gr) %>% as.character() %>% unique() %>%  gtools::mixedsort()
    seqlevels(region.gr.shift) <- seqnames(region.gr.shift) %>% as.character() %>% unique() %>%  gtools::mixedsort()
    chunk.l <- scanBam(file = bam, param = ScanBamParam(tag = tag,
                                                        what = bamFanc::ALL.FIELDS,
                                                        which = region.gr.shift + padding))
    chunk.l <- bamFanc::bam.chunk.to.df.2(chunks = chunk.l, rbind = F,
                                          add.bits = c("isSecondaryAlignment", "isFirstMateRead"),
                                        add.min.mapq = T, add.m.score = T)
    ## rename chunk.l
    rename.vec <- mcols(region.gr.shift)[, name.col]
    names(rename.vec) <- utilsFanc::gr.get.loci(region.gr.shift + padding)
    if (!identical(names(chunk.l), names(rename.vec)))
      stop("!identical(names(chunk.l), names(rename.vec))")
    names(chunk.l) <- rename.vec[names(chunk.l)]
    # chunk.l <- bamFanc::bam.df.subset.by.region(bam.df = chunk, gr = shift(region.gr, -shift),
    #                                             return.list = T, gr.name.col = name.col,
    #                                             rescue.qname = F)
    chunk.l.f <- bamFanc::nkc.composite.filter(chunk.l = chunk.l, m.score.cutoff = m.score.cutoff,
                                  out.dir = paste0(report.dir, "/", sample,"/"),
                                  threads = threads.region, header.file = header.file,
                                  tags = tag, shift = shift)

    chunk.l.f <- lapply(chunk.l.f, function(x) return(x$pass))
    cutsite.l <- bamFanc::bam.get.cut.site.ez(chunk.l = chunk.l.f, overall.shift = shift,
                                              left.shift = 0, right.shift = 0,
                                              gr.region.name.map = region.gr + padding,
                                              name.col = name.col)
    return(cutsite.l)
  }, mc.cores = threads.bam)
  names(cutsite) <- sample.names
  stat <- lapply(cutsite, function(x) {
    n.cut <- sapply(x, function(y) return(length(y$gr)))
    return(n.cut)
  }) %>% Reduce(cbind, .) %>% as.data.frame() %>% `colnames<-`(names(cutsite)) %>%
    mutate(., region = rownames(.)) %>% `rownames<-`(NULL)
  write.table(stat, paste0(report.dir, "/stat.tsv"), quote = F,sep = "\t", row.names = F, col.names = T)

  if (!is.null(ctrl.bed.files)) {
    names(ctrl.bed.files) <- sample.names
    cutsite.ctrl <- mclapply(sample.names, function(sample) {
      bed <- ctrl.bed.files[sample]
      bed.l <- (region.gr + padding) %>%
        split(., f = factor(mcols(.)[, name.col], levels = unique(mcols(.)[, name.col]))) %>%
        mclapply(function(gr) {
          bed.gr <- rtracklayer::import.bed(bed, which = gr)
        }, mc.cores = threads.region, mc.cleanup = T)
      cutsite.l <- bamFanc::bam.get.cut.site.ez(chunk.l = bed.l, is.open.bed = T,
                                                overall.shift = 0,
                                                # left.shift = 0, right.shift = 0,
                                                gr.region.name.map = region.gr,
                                                name.col = name.col)
      return(cutsite.l)
    }, mc.cores = threads.bam, mc.cleanup = T)
    names(cutsite.ctrl) <- sample.names
  }
  stat.ctrl <- lapply(cutsite.ctrl, function(x) {
    n.cut <- sapply(x, function(y) return(length(y$gr)))
    return(n.cut)
  }) %>% Reduce(cbind, .) %>% as.data.frame() %>% `colnames<-`(names(cutsite)) %>%
    mutate(., region = rownames(.)) %>% `rownames<-`(NULL)

  cutsite.join <- mapply(function(fg, ctrl) {
    union <- mapply(function(x, y) {
      common.qnames <- intersect(x$df$qname, y$df$qname)
      x$df <- rbind(x$df %>% mutate(type = "fg"), y$df %>% mutate(type = "ctrl") %>%
                    filter(!qname %in% common.qnames) %>% dplyr::select(chr, start, end, qname, type))
      x$gr <- NULL
      x$df.pe <- NULL
      return(x)
    }, fg, ctrl, SIMPLIFY = F)
    return(union)
  }, cutsite, cutsite.ctrl, SIMPLIFY = F)

  stat.join <- lapply(cutsite.join, function(x) {
    n.cut <- sapply(x, function(y) return(nrow(y$df)))
    return(n.cut)
  }) %>% Reduce(cbind, .) %>% as.data.frame() %>% `colnames<-`(names(cutsite)) %>%
    mutate(., region = rownames(.)) %>% `rownames<-`(NULL)
  res <- list(fg = list(cutsite = cutsite, stat = stat),
              ctrl = list(cutsite = cutsite.ctrl, stat = stat.ctrl),
              join = list(cutsite = cutsite.join, stat = stat.join))
  return(res)
}
