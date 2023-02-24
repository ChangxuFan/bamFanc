bam.shift.file <- function(bam, bam.out = NULL, shift.chr = NULL, shift = 0,
                           samtools = "/opt/apps/samtools/1.9/bin/samtools") {
  if (shift == 0)
    return(bam)
  if (is.null(bam.out)) {
    bam.out <- bam
  }
  tmp <- tempfile()
  matching <- "$1 !~ /^@/"
  if (is.null(shift.chr)) {
    matching <- paste0(matching, " && $3 == ", shift.chr)
  }
  cmd <- paste0("samtools view -h ", bam, " | awk -F '\\t' 'BEGIN {OFS = FS} ",
                matching,
                "{$4 = $4 + ", shift, "; print $0}", " ' ")
}