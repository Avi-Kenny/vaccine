#' Run diagnostics
#'
#' @description TO DO
#' @param obj An object of class \code{vaccine_est} returned by est_ce
#' @return A plot of model diagnostics
#' @examples
#' print("to do")
#' @export
diagnostics <- function(obj) {

  if (!methods::is(obj,"vaccine_est")) {
    stop(paste0("`obj` must be an object of class 'vaccine_est' returned by es",
                "t_ce()."))
  }

  p1 <- ggplot2::ggplot(obj$extras$r_Mn, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="r_Mn")
  p2 <- ggplot2::ggplot(obj$extras$deriv_r_Mn, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="deriv_r_Mn")
  p3 <- ggplot2::ggplot(obj$extras$Gamma_os_n, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="Gamma_os_n")
  p4 <- ggplot2::ggplot(obj$extras$f_s_n, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="f_s_n")

  p5 <- ggplot2::ggplot(
    obj$extras$Q_n,
    ggplot2::aes(x=t, y=est, group=factor(interaction(ind,s)), color=factor(s))
  ) +
    ggplot2::geom_line(alpha=0.5) +
    ggplot2::facet_wrap(~s) +
    ggplot2::labs(y="", title="Q_n", color="s") +
    ggplot2::theme(legend.position="none")

  p6 <- ggplot2::ggplot(
    obj$extras$Qc_n,
    ggplot2::aes(x=t, y=est, group=factor(interaction(ind,s)), color=factor(s))
  ) +
    ggplot2::geom_line(alpha=0.5) +
    ggplot2::facet_wrap(~s) +
    ggplot2::labs(y="", title="Qc_n", color="s") +
    ggplot2::theme(legend.position="none")

  plot <- ggpubr::ggarrange(plotlist=list(p5,p1,p2,p6,p3,p4), ncol=3, nrow=2)

  return(plot)

}
