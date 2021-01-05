usethis::use_package('dplyr')
usethis::use_package('ggplot2')

wideToLongCIandEst <- function(ci.df, est.df, truth) {
  n.contrast <- nrow(ci.df$polyhedral.symmetric)
  ci.list    <- lapply(ci.df, unlist)
  ci.type    <- names(ci.list)
  ci.df.long <- data.frame('Type' = NULL, 'Lower' = NULL, 'Upper' = NULL)
  for (i in 1:length(ci.type)) {
    ci.df.long <- suppressWarnings(dplyr::bind_rows(ci.df.long,
                                                    data.frame('Index' = seq(1, n.contrast, 1),
                                                               'Type'  = rep(ci.type[i], n.contrast),
                                                               'Lower' = ci.list[[i]][ ,1],
                                                               'Upper' = ci.list[[i]][ ,2])))
  }
  est.df['Truth'] <- truth
  if (any(colnames(est.df) %in% 'moment')) {
    est.df <- dplyr::rename(est.df, 'polyhedral.symmetric' = 'moment')
  }
  est.type    <- names(est.df)
  est.df.long <- data.frame('Type' = NULL, 'Est' = NULL)
  for (i in 1:length(est.type)) {
    est.df.long <- suppressWarnings(dplyr::bind_rows(est.df.long,
                                                     data.frame('Index'    = seq(1, n.contrast, 1),
                                                                'Type'     = rep(est.type[i], n.contrast),
                                                                'Estimate' = est.df[[i]])))
  }
  plot.df <- dplyr::full_join(est.df.long, ci.df.long)
  return(plot.df)
}


#' @title Plotting results of PSAT procedures
#'
#' @description \code{plotFunc} is used for basic plotting of the CI and estimates resulting from the PSAT procedure.
#' @param ci.df data.frame of confidence intervals chosen to plot
#' @param est.df data.frame of estimates chosen to plot, `moment` estimates name is changed to `polyhedral.symmetric`
#' @param base.line Draw horiziontal line
#' @details Offers some basic plotting of the
#' @examples
#' \dontrun{
#'
#' y         <- c(15, rnorm(p - 1))
#' cov.mat <- diag(p)
#' psat.res <- PSATQuadratic(y, cov.mat) ### Runs on default values
#' plotFunc(psat.res$CI %>%
# '         select(global, switch.symmetric, naive, polyhedral.symmetric), psat.res$Point.estimation))
#'
#' }
#'
#' @export
plotFunc <- function(ci.df, est.df, truth = NULL, base.line = NULL) {
  plot.df <- wideToLongCIandEst(ci.df, est.df, truth)
  ci.plot <- suppressWarnings(ggplot2::ggplot(plot.df) +
    ggplot2::geom_linerange(aes(ymin = Lower,
                       ymax = Upper,
                       x = Index,
                       color = Type),
                   position = position_dodge(0.5)) +
    ggplot2::geom_point(aes(x = Index,
                   y = Estimate,
                   shape = Type,
                   color = Type),
               position = position_dodge(0.5)))
    if (!is.null(base.line)) {
      ci.plot <- ci.plot + ggplot2::geom_segment(aes(y = base.line, yend = base.line, x = 0, xend = max(plot.df$Index)))
    }
  return(ci.plot)
}
