# rxy_jackknife.R — Rxy relative load with equal-block weighted jackknife
#
# Rxy statistic: Do et al. 2015, Nat Genet 47:126-131
# Weighted block jackknife SE: Busing et al. 1999, Stat Comput 9:3-8
# Block generation and Busing functions from Stuart 2026, Mol Ecol Resour 
# code adapted from: https://github.com/OliverPStuart/2025_Jackknife_Review.git
#
#
# Input:
#   site_class.txt       CHROM, POS, IMPACT, EFFECT (from SnpSift extractFields)
#   af_Ire_nocurl1.txt   Per-population derived (ALT) allele frequencies
#   af_Brit.txt          Format: CHROM <tab> POS <tab> AF
#   af_Swe.txt
#   lengths.txt          Scaffold name <tab> length in bp
#
# Output:
#   rxy_results_equalblocks.txt

setwd(".")
library(data.table)

# Stuart 2026 / Busing 1999 functions,
weighted_jackknife_se <- function(x, m){
  n <- sum(m); theta_hat <- weighted.mean(x, m)
  sqrt(mean(((n - m)/m) * (theta_hat - x)^2))
}
weighted_jackknife_mean <- function(x, m, theta_hat){
  g <- length(x); n <- sum(m)
  g*theta_hat - sum(((n - m)*x)/n)
}

# data
lengths <- fread("lengths.txt", header=FALSE); setnames(lengths, c("Scaffold","Length"))
site <- fread("site_class.txt"); setnames(site, c("chrom","pos","impact","effect"))
site[, class := fifelse(impact=="MODIFIER","NEUTRAL",
                 fifelse(impact=="MODERATE","MODERATE",
                 fifelse(impact=="HIGH","HIGH", NA_character_)))]
site <- site[!is.na(class), .(chrom,pos,class)]
readaf <- function(p){ d <- fread(paste0("af_",p,".txt"), na.strings="."); setnames(d,c("chrom","pos","af")); d[!is.na(af)] }
af <- list(Ire=readaf("Ire_nocurl1"), Brit=readaf("Brit"), Swe=readaf("Swe"))

# Stuart 2026 block generation (physical blocks with redistribution)
n_blocks <- 101
total_len <- sum(lengths$Length)
block_size <- total_len / n_blocks
blocks <- data.frame()
for (i in seq_len(nrow(lengths))) {
  chr <- lengths$Scaffold[i]; chr_len <- lengths$Length[i]
  n_chr <- max(1, round(chr_len / block_size)); chr_block_size <- chr_len / n_chr
  edges <- round(seq(1, chr_len + 1, length.out = n_chr + 1))
  starts <- edges[-length(edges)]; ends <- edges[-1] - 1
  df <- data.frame(Scaffold = chr, Start = starts, End = ends)
  k <- nrow(df)
  if (k >= 2) {
    last_size <- df$End[k] - df$Start[k] + 1L
    frac <- last_size / chr_block_size
    if (last_size > 0 && frac <= 0.75) {
      if (frac <= 0.25) m <- 2 else m <- 3
      m <- min(m, k - 1L)
      sizes <- df$End - df$Start + 1L
      add_base <- last_size %/% m; rem <- last_size - add_base * m
      idx <- seq.int(k - m, k - 1)
      sizes[idx] <- sizes[idx] + add_base
      if (rem > 0L) sizes[tail(idx, rem)] <- sizes[tail(idx, rem)] + 1L
      sizes <- sizes[-k]
      new_starts <- cumsum(c(1, head(sizes, -1))); new_ends <- cumsum(sizes)
      df <- data.frame(Scaffold = chr, Start = new_starts, End = new_ends)
    }
  }
  blocks <- rbind(blocks, df)
}
blocks$Block <- 1:nrow(blocks)
setDT(blocks)

## --- SNP -> block assignment: foverlaps (identical result to Stuart's per-row filter) ---
assign_blocks <- function(d){
  d2 <- copy(d); d2[, `:=`(start=pos, end=pos)]
  setkey(blocks, Scaffold, Start, End)
  ov <- foverlaps(d2, blocks,
                  by.x=c("chrom","start","end"),
                  by.y=c("Scaffold","Start","End"),
                  type="within", nomatch=NA)
  ov$Block
}

# Rxy per pair, standardised to NEUTRAL, jackknife over Stuart's blocks ---
rxy_pair <- function(A, B){
  d <- merge(merge(site, af[[A]], by=c("chrom","pos")),
             af[[B]], by=c("chrom","pos"), suffixes=c(".A",".B"))
  d[, LXY := af.A*(1-af.B)][, LYX := af.B*(1-af.A)]
  d <- d[is.finite(LXY) & is.finite(LYX)]
  d[, Block := assign_blocks(d)]
  d <- d[!is.na(Block)]

  rxy_cls <- function(dt){
    r <- dt[, .(RXY=sum(LXY)/sum(LYX)), by=class]
    r[, RXY_st := RXY/r[class=="NEUTRAL", RXY]][]
  }
  whole <- rxy_cls(d)

  blks <- sort(unique(d$Block))
  m <- sapply(blks, function(b) d[Block==b, .N])

  bs  <- d[, .(sLXY=sum(LXY), sLYX=sum(LYX)), by=.(Block, class)]
  tot <- bs[, .(TLXY=sum(sLXY), TLYX=sum(sLYX)), by=class]
  ps_for <- function(b, cl){
    left <- merge(tot, bs[Block==b, .(class,sLXY,sLYX)], by="class", all.x=TRUE)
    left[is.na(sLXY),sLXY:=0][is.na(sLYX),sLYX:=0]
    left[, RXY:=(TLXY-sLXY)/(TLYX-sLYX)]
    neut <- left[class=="NEUTRAL", RXY]
    left[class==cl, RXY/neut]
  }

  rbindlist(lapply(c("MODERATE","HIGH"), function(cl){
    theta <- whole[class==cl, RXY_st]
    ps <- sapply(blks, ps_for, cl=cl)
    mean_j <- weighted_jackknife_mean(ps, m, theta)
    se_j   <- weighted_jackknife_se(ps, m)
    g <- length(blks); t <- (mean_j - 1)/se_j
    data.table(pair=paste(A,B,sep="_"), class=cl,
               Rxy=round(mean_j,3), SE=round(se_j,3),
               t=round(t,2), p=round(2*(1-pt(abs(t), g-1)),4), n_blocks=g)
  }))
}

res <- rbindlist(lapply(list(c("Ire","Brit"),c("Ire","Swe"),c("Brit","Swe")),
                        function(p) rxy_pair(p[1], p[2])))
print(res)
fwrite(res, "rxy_results_equalblocks.txt", sep="\t")
