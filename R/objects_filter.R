bam.filter.chunk <- function(chunk, eof, bof, hinge.pre, filter.fun, filter.fun.params = NULL,
                             outfile, junkfile = NULL, stat.file = NULL, header.file, 
                             add.index = T, threads = 6,
                             samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools") {
  if (is.null(junkfile))
    junkfile <- utilsFanc::insert.name(name = outfile, insert = "junk",
                                       ext = ".bam", trim.dir = F)
  final.files  <- c(outfile = outfile, junkfile = junkfile)

  final.files.w.mate <- utilsFanc::insert.name(final.files, insert = "withMate",
                                               ext = ".bam", trim.dir = F)

  tempt.files <- paste0(final.files, ".sam")

  if (eof == T) {
    lapply(1:2, function(i) {
      cmd <- paste0(samtools," view -b -h -@ ", threads, " ", tempt.files[i], " > ", final.files.w.mate[i])
      utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)
      if (add.index == T)
        system(paste0(samtools, " index ", final.files.w.mate[i]))

      if (names(final.files)[i] != "junkfile") {
        remove.mate(bam = final.files.w.mate[i], out.bam = final.files[i], thread = threads)
        if (add.index == T)
          system(paste0(samtools, " index ", final.files[i]))
      }
      return(NULL)
    })
    if (is.null(stat.file))
      stat.file <- sub(".bam$",".stat",outfile)

    files.to.count <- c(final.files.w.mate, outfile)
    counts <- mclapply(files.to.count, function(file) {
      cmd <- paste0(samtools, " view -c -@ 6 ", file)
      n <- utilsFanc::cmd.exec.fanc(cmd = cmd, intern = T, run = T)
      return(n)
    }, mc.cores = length(files.to.count), mc.cleanup = T) %>% unlist()
    stat <- data.frame(files = basename(files.to.count),
                       nreads = counts)
    write.table(stat, stat.file, col.names = T, row.names = F, quote = F, sep = "\t")
    return(NULL)
  }
  # append.mode <- T
  if (bof == T) {
    lapply(1:2, function(i) {
      system(paste0("cp ", header.file, " ", tempt.files[i]))
    })
  }
  # first perform some basic filtering
  index <- do.call(what = filter.fun, args =c(list(chunk = chunk), filter.fun.params))
  if (!is.logical(index)) {
    stop("filter.fun must return a logical vector")
  }
  df.o <- chunk[index, ] %>% scFanc::factor2character.fanc()
  df.j <- chunk[!index, ] %>% scFanc::factor2character.fanc()
  dfs <- list(df.o = df.o, df.j = df.j)
  lapply(1:2, function(i) {
    df <- dfs[[i]]
    # browser()
    if (nrow(df) >= 1 ) {
      df <- df %>% tag.formatting() %>% rnext.formatting()
      write.table(df, tempt.files[i], append = T, quote = F, col.names = F,
                  row.names = F, sep = "\t")
    }
  })
  return(NULL)
}

bam.filter.no.mismatch <- function(chunk, tlen.cutoff = NULL) {
  bIndex.mismatch <- !grepl("[ATCG]", chunk$tag.MD)
  bIndex.indel <- !grepl("[ID]", chunk$cigar)
  bIndex <- bIndex.mismatch & bIndex.indel
  if (!is.null(tlen.cutoff)) {
    bIndex.tlen <- abs(chunk$isize) < tlen.cutoff
    bIndex <- bIndex & bIndex.tlen
  }
  return(bIndex)
}

bam.filter.softclip <- function(chunk, general.cutoff, tss.cutoff = 3, threads = 6) {
  bIndex <- chunk %>% split(f = 1:nrow(chunk)) %>%
    mclapply(function(x) {
      # print(x$qname)
      s.length <- explodeCigarOpLengths(x$cigar, "S") %>% unlist()
      # print(paste0("s.length: ",s.length))
      if (sum(s.length > general.cutoff) > 0)
        return(F)

      flags <- bamFlagAsBitMatrix(x$flag)
      if (flags[1, "isFirstMateRead"] == 0) {
        # print("not first mate")
        return(T)
      }

      # print(x$cigar)
      cigar.op <- explodeCigarOps(x$cigar) %>% unlist()
      cigar.l <- explodeCigarOpLengths(x$cigar) %>% unlist()

      if (flags[1, "isMinusStrand"] == 0) {
        # pos strand, look at the first
        pointer <- 1
        nuc <- "G"
      } else {
        pointer <- length(cigar.op)
        nuc <- "C"
      }
      # print(paste0("pointer: ", pointer, "; nuc: ",nuc))
      # print(x$seq)
      # if there is no softclipping at the pointer side, keep the read:
      if (cigar.op[pointer] != "S")
        return(T)
      if (cigar.l[pointer] > tss.cutoff)
        return(F)

      if (nuc == "G") {
        clipped <- substr(x$seq, 1, cigar.l[pointer])
        # print("clipped: ")
        # print(clipped)
        if (grepl("[^G]", clipped)) {
          # print("F")
          return(F)
        }
        else {
          # print("T")
          return(T)
        }
      } else {
        clipped <- substr(x$seq, nchar(x$seq) - cigar.l[pointer] + 1, nchar(x$seq))
        # print("clipped: ")
        # print(clipped)
        if (grepl("[^C]", clipped)) {
          # print("F")
          return(F)
        }
        else {
          # print("T")
          return(T)
        }
      }
    }, mc.cores = threads, mc.cleanup = T) %>% unlist()
  names(bIndex) <- NULL
  return(bIndex)
}


remove.mate <- function(bam, out.bam = NULL, debug = F,
                        samtools = "/bar/cfan/anaconda2/envs/jupyter/bin/samtools",
                        sambamba = "/bar/cfan/software/sambamba/sambamba",
                        thread = 6) {
  if (is.null(out.bam))
    out.bam <- utilsFanc::insert.name(name = bam, insert = "removeMate",
                                      ext = ".bam", trim.dir = F)
  bam.nsort <- utilsFanc::insert.name(name = bam, insert = "nsort",
                           ext = ".bam", trim.dir = F)
  cmd <- paste0(samtools, " sort -n ", " -@ ", thread,
                " -o ", bam.nsort, " ", bam)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)
  bam.nsort.fixmate <- utilsFanc::insert.name(name = bam.nsort, insert = "fixmate",
                                   ext = ".bam", trim.dir = F)
  cmd <- paste0(samtools, " fixmate ", bam.nsort, " ",
                bam.nsort.fixmate)
  utilsFanc::cmd.exec.fanc(cmd = cmd, intern = F, run = T)

  cmd <- paste0(samtools, " sort ", " -@ ", thread, " ",
                bam.nsort.fixmate, " | ", sambamba, " view -f bam -F ",
                "\"", "paired", "\"", " -t ", thread, " ", "/dev/stdin",
                " > ", out.bam,
                " && ", samtools, " index ", out.bam)
  print(cmd)
  system(cmd)
  if (debug == F) {
    system(paste0("rm ", bam.nsort, " ", bam.nsort.fixmate))
  }
  if (!file.exists(out.bam))
    stop(paste0(out.bam, " was not successfully generated"))
  return(out.bam)
}

bam.filter.mqa <- function(bam, out.bam = NULL, valid.reads) {
  reads <- readLines(valid.reads) %>% unique()
  filter <- function(x) x$qname %in% reads
  filters <- list(name.filter = filter)
  if (is.null(out.bam))
    out.bam <- sub(".bam", "_mqa.bam", bam)
  Rsamtools::filterBam(file = bam, destination = out.bam, filter = FilterRules(filters))
  return(out.bam)
}