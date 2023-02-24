one.mate.multimap <- function(chunk.df, strict = F, threads = 12) {
  if (is.null(chunk.df$bFirst))
    chunk.df$bFirst <- bamFlagAsBitMatrix(chunk.df$flag)[, "isFirstMateRead"]
  if (strict == T) {
    bOmm <- chunk.df %>% split(f = chunk.df$qname) %>% 
      mclapply(function(df) {
        if (any(is.na(df$pos)))
          return(F)
        if (min(table(df$bFirst)) == 1) {
          if (max(table(df$bFirst)) == 1)
            return(F)
          else
            return(T)
        } else {
          AS <- df$tag.AS %>% split(df$bFirst) %>% 
            sapply(function(x) {
              return(sum(x == max(x)))
            })
          if (min(AS) == 1 && max(AS) > 1)
            return(T)
          else
            return(F)
        }
      }, mc.cores = threads, mc.cleanup = T)
    qname.omm <- names(bOmm)[bOmm == T]
    qname.not.omm <- names(bOmm)[bOmm == F]
    omm.rate <- round(length(qname.omm)/length(qname.not.omm), digit = 3)
    filtered.df <- chunk.df %>% filter(qname %in% qname.omm)
    res <- list(qname.omm = qname.omm, omm.rate = omm.rate, filtered.df = filtered.df)
  } else {
    occur.df <- chunk.df %>% group_by(qname) %>% 
      summarise(n.mate1 = sum(bFirst == 1), n.mate2 = sum(bFirst == 0))
    occur.df$f.mate1 <- 1
    occur.df$f.mate1[occur.df$n.mate1 != 1] <- -1
    
    occur.df$f.mate2 <- 1
    occur.df$f.mate2[occur.df$n.mate2 != 1] <- -1
    
    occur.df <- occur.df %>% mutate(x = f.mate1 * f.mate2)
    qname.omm <- occur.df$qname[occur.df$x == -1]
    
    filtered.df <- chunk.df %>% filter(qname %in% qname.omm)
    res <- list(qname.omm = qname.omm, filtered.df = filtered.df, occur.df = occur.df)
  }
  
  return(res)
}