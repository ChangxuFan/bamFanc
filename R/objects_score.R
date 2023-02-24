score.distro <- function(bams, sample.names, out.dir, threads = NULL, pub = T) {
  if (is.null(threads))
    threads <- length(bams)
  dir.create(out.dir, showWarnings = F, recursive = T)
  res <- utilsFanc::safelapply(seq_along(sample.names), function(i) {
    sample <- sample.names[i]
    bam <- bams[i]
    bam.df <- Rsamtools::scanBam(file = bam, param = ScanBamParam(what = c("qwidth", "qname"), tag = c("AS", "XS"))) %>% 
      bam.chunk.to.df.2(rbind = T)
    
    bam.df <- bam.df %>% mutate(score = tag.AS/qwidth, sample = sample)
    stat <- data.frame(sample = sample, n.aln = nrow(bam.df), 
                       mean = mean(bam.df$score), sd = sd(bam.df$score),
                       q10 = quantile(bam.df$score, 0.1))
    stat$min_sd <- stat$mean - stat$sd
    q <- quantile(x = bam.df$score, 0.05*(1:20))
    q.df <- data.frame(quantile = names(q), value = q)
    write.table(stat, paste0(out.dir, "/", sample, "_stat.tsv"), sep = "\t", quote = F, row.names = F,
                col.names = T)
    write.table(q.df, paste0(out.dir, "/", sample, "_quantile.tsv"), sep = "\t", quote = F, 
                row.names = F, col.names = T)
    if (pub) {
      subset.df <- bam.df[sort(sample(1:nrow(bam.df), 100000)),]
      p <- ggplot(subset.df, aes(x = score)) +
        geom_density(adjust = 0.5, color = "orangered", fill = alpha("orangered", 0.3), size = 0.3) +
        geom_vline(xintercept = stat$q10, color = "blue", linetype = "dotted", size = 0.4)
      p <- p %>% utilsFanc::theme.fc.1(italic.x = F) + theme(aspect.ratio=1)
      ggsave(paste0(out.dir, "/", sample, "_distro.pdf"), p, device = cairo_pdf, 
             width = 1.5, height = 1.5, dpi = 300, units = "in")
      
    } else {
      p <- ggplot(bam.df, aes(x = score)) +
        geom_density() +
        geom_vline(xintercept = c(stat$mean), color = "blue", linetype="dotted") +
        geom_vline(xintercept = c(stat$min_sd), color = "grey20", linetype = "dotted") +
        theme_bw() +
        theme(aspect.ratio=1) + 
        ggsave(paste0(out.dir, "/", sample, "_distro.png"), width = 7, height = 5, dpi = 100, units = "in")
      res <- list(stat = stat, quantile = q.df, p = p)
      return(res)
    }
  }, threads = threads)
  names(res) <- sample.names
  return(res)
}