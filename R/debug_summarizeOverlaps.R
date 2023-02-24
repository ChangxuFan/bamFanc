# samples.df <- data.frame(bamfile = c("~/4dn/nk/oshea_2020/sth/raw/histone_refilter/K4me3_NK_WT_exvivo_rep1_chip_GSM4314396_oshea2021_SRR11087490_fastp_mkdup_filter_mapq_1.bam",
#                                      "~/4dn/nk/oshea_2020/sth/raw/histone_refilter/K4m3_NK_WT_exvivo_rep2_chip_GSM4314407_oshea2032_SRR11087501_fastp_mkdup_filter_mapq_1.bam"),
#                          sample = c("rep1", "rep2"))
# 
# gr.mini <- utilsFanc::loci.2.df(loci.vec = c("chr1:150,394,298-150,394,391",
#                                              "chr1:150,394,298-150,394,482"),
#                                 from.igv = T, remove.loci.col = T, return.gr = T)
# 
# strand(gr.mini) <- "+"
# 
# bam1 <- BamFile(samples.df$bamfile[1], index = paste0(samples.df$bamfile[1]))
# a <- Rsamtools::scanBam(file = bam1,
#                         param = ScanBamParam(which = gr.mini, what = bamFanc::ALL.FIELDS))
# ## takes a second to finish
# tt <- GenomicAlignments::summarizeOverlaps(features = gr.mini, reads = bam1, 
#                                            mode = "Union", inter.feature = F,
#                                            singleEnd = bSingleEnd, fragments = F,
#                                            ignore.strand = F)
# ## takes forever. It seems that summarizeOverlaps will read the entire bam file
# b <- readGAlignments(samples.df$bamfile[1], use.names=FALSE, param=ScanBamParam(which = gr.mini),
#                      with.which_label=FALSE)
# ## very fast
# d <- GenomicAlignments::summarizeOverlaps(features = gr.mini, reads = list(miao = b, wang = b), 
#                                           mode = "Union", inter.feature = F,
#                                           singleEnd = bSingleEnd, fragments = F,
#                                           ignore.strand = F)
# ## doesn't work
# 
# bam2 <- BamViews(samples.df$bamfile, bamRanges = gr.mini)
# b2 <- readGAlignments(bam2, use.names=FALSE, param=ScanBamParam(which = gr.mini),
#                       with.which_label=FALSE)
# b3 <- do.call(GAlignmentsList, b2)
# 
# e <- GenomicAlignments::summarizeOverlaps(features = gr.mini, reads = b3[[1]],
#                                           mode = "Union", inter.feature = F,
#                                           singleEnd = bSingleEnd, fragments = F,
#                                           ignore.strand = T)
# 
# ## each element in b3 is treated as one interval with multiple segments