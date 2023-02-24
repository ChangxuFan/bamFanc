bam2fq <- function(bams, regions, root.names = NULL, out.dir, fq.suffix = ".fastq.gz", add.R = "T",
                   threads.bam = 1, threads.region = 1, threads.samtools = 1,
                   samtools = "/opt/apps/samtools/1.9/bin/samtools") {
  # regions: a list of grs. or just 1 gr.
  if (is.null(root.names)) {
    root.names <- basename(bams) %>% sub(".bam$", "", .)
  }
  if (!is.list(regions) || is.null(names(regions))) {
    stop("regions must be a named list of GRanges!")
  }
  fastq.list <- utilsFanc::safelapply(seq_along(bams), function(i) {
    bam <- bams[i]
    root.name <- root.names[i]

    fastq.by.region <- utilsFanc::safelapply(names(regions), function(region.name) {
      gr <- regions[[region.name]]
      if (! "GRanges" %in% class(gr)) {
        stop(! "GRanges" %in% class(gr))
      }
      loci <- utilsFanc::gr.get.loci(gr = gr) %>% paste0(collapse = " ")
      if (add.R) {
        add.R <- "_R"
      } else {
        add.R <- "_"
      }
      fastqs <- paste0(out.dir, "/", root.name, "_", region.name, add.R, 1:2, fq.suffix)
      system(paste0("mkdir -p ", out.dir))
      cmd <- paste0(samtools, " view -hb ", bam, " ", loci, " | ",
                    samtools, " sort -n -m 4G -@ ", threads.samtools, " - | ",
                    samtools, " fastq -n -1 ", fastqs[1], " -2 ", fastqs[2],
                    " -0 /dev/null -s /dev/null -")
      print(cmd)
      system(cmd)
      return(fastqs)
    }, threads = threads.region)
    names(fastq.by.region) <- names(regions)
    return(fastq.by.region)
  }, threads = threads.bam)
  names(fastq.list) <- root.names
  return(fastq.list)
}