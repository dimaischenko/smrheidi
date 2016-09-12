#' Function to calculate SMR statistic
#' 
#' @param gw.pv vector with GWAS p-values
#' @param eq.pv vector with eQTL p-values
#' 
#' @return SMR p-values
smr_test <- function(gw.pv, eq.pv) {
  # calc z-value for GWAS
  gw.z <- qnorm(p = gw.pv / 2, lower.tail = F)
  # calc z-value for eQTL
  eq.z <- qnorm(p = eq.pv / 2, lower.tail = F)
  
  # calc SMR statistic
  t_smr <- (gw.z^2 * eq.z^2) / (gw.z^2 + eq.z^2)
  
  # return p-value
  return(pchisq(t_smr, df = 1, lower.tail = F))
}

#' Function to calculate HEIDI statistic
#'
#' @param gw.b vector with beta for GWAS snps
#' @param gw.bse vector with beta standard errors for GWAS snps
#' @param eq.b vector with beta for eQTL snps
#' @param eq.bse vector with beta standard errors for eQTL snps
#' @param r vector with LD correlation between top snp and all snps in region
#' @param top numeric index of top snp in vectors
#' 
#' @return HEIDI p-value
heidi_test <- function(gw.b, gw.bse, eq.b, eq.bse, r, top) {
  # lib for Satterthwaite's approximation
  require(survey)
  
  # gwas z-values
  gw.z <- gw.b / gw.bse
  eq.z <- eq.b / eq.bse
  
  # calc b_xy
  b <- gw.b / eq.b
  # calc d
  d <- b - b[top]
  
  # calc b standard error
  b.se <- sqrt(
    gw.b^2/eq.b^2 * (gw.bse^2/gw.b^2 + eq.bse^2/eq.b^2)
  )
  
  # covariance with top snp
  b.cov.top <- r / (eq.b * eq.b[top]) * gw.bse * gw.bse[top] + 
    b * b[top] * (r/(eq.z * eq.z[top]) - 1/(eq.z^2 * eq.z[top]^2))
  
  # variance of d
  d.se <- b.se^2 + b.se[top]^2 - 2 * b.cov.top
  
  # calc z-values for {d}
  d.z <- d[-top] / sqrt(d.se[-top])
  
  # TODO(dima): check why b.se can be negative, is it ok?
  # now we exclude this snps from analysis
  d.z <- d.z[!is.na(d.z)]
  
  # calc HEIDI statistic
  t_heidi <- sum(d.z^2)
  
  # calc p-value
  return(pchisqsum(t_heidi, df = rep(1, length(d.z)),
                   a = rep(1, length(d.z)),
                   lower.tail = F,  method = "satterthwaite"))
}
