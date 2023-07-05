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

  p1 <- ggplot(obj$extras$r_Mn, aes(x=s, y=est)) +
    geom_line(color="#56B4E9") + geom_point(color="#56B4E9", size=1) +
    labs(y="", x="log(S)", title="r_Mn")
  p2 <- ggplot(obj$extras$deriv_r_Mn, aes(x=s, y=est)) +
    geom_line(color="#56B4E9") + geom_point(color="#56B4E9", size=1) +
    labs(y="", x="log(S)", title="deriv_r_Mn")
  p3 <- ggplot(obj$extras$Gamma_os_n, aes(x=s, y=est)) +
    geom_line(color="#56B4E9") + geom_point(color="#56B4E9", size=1) +
    labs(y="", x="log(S)", title="Gamma_os_n")
  p4 <- ggplot(obj$extras$f_s_n, aes(x=s, y=est)) +
    geom_line(color="#56B4E9") + geom_point(color="#56B4E9", size=1) +
    labs(y="", x="log(S)", title="f_s_n")

  p5 <- ggplot(
    obj$extras$Q_n,
    aes(x=t, y=est, group=factor(interaction(ind,s)), color=factor(s))
  ) +
    geom_line(alpha=0.5) +
    facet_wrap(~s) +
    labs(y="", title="Q_n", color="s") +
    theme(legend.position="none")

  p6 <- ggplot(
    obj$extras$Qc_n,
    aes(x=t, y=est, group=factor(interaction(ind,s)), color=factor(s))
  ) +
    geom_line(alpha=0.5) +
    facet_wrap(~s) +
    labs(y="", title="Qc_n", color="s") +
    theme(legend.position="none")

  plot <- ggpubr::ggarrange(plotlist=list(p5,p1,p2,p6,p3,p4), ncol=3, nrow=2)

  return(plot)

}
