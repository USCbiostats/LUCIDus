#' @title Visualize LUCID model through a Sankey diagram
#' @description In the Sankey diagram, each node either represents a variable (exposure,
#' omics or outcome) or a latent cluster. Each line represents an association. The
#' color of the node represents variable type, either exposure, omics or outcome.
#' The width of the line represents the effect size of a certain association; the
#' color of the line represents the direction of a certain association. 
#' 
#' @param x A LUCID model fitted by \code{\link{est_lucid}}
#' @param G_color Color of node for exposure
#' @param X_color Color of node for latent cluster
#' @param Z_color Color of node for omics data
#' @param Y_color Color of node for outcome
#' @param pos_link_color Color of link corresponds to positive association
#' @param neg_link_color Color of link corresponds to negative association
#' @param fontsize Font size for annotation
#' 
#' @return A DAG graph created by \code{\link{sankeyNetwork}}
#' 
#' @import networkD3
#' @importFrom jsonlite toJSON
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # prepare data
#' G <- sim_data$G
#' Z <- sim_data$Z
#' Y_normal <- sim_data$Y_normal
#' Y_binary <- sim_data$Y_binary
#' cov <- sim_data$Covariate
#' 
#' # plot lucid model
#' fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", 
#' K = 2, seed = 1008)
#' plot_lucid(fit1)
#' 
#' # change node color
#' plot_lucid(fit1, G_color = "yellow")
#' plot_lucid(fit1, Z_color = "red")
#' 
#' # change link color
#' plot_lucid(fit1, pos_link_color = "red", neg_link_color = "green")
#' }
plot_lucid <- function(x,
                       G_color = "dimgray",
                       X_color = "#eb8c30",
                       Z_color = "#2fa4da",
                       Y_color =  "#afa58e",
                       pos_link_color = "#67928b", 
                       neg_link_color = "#d1e5eb",
                       fontsize = 7
                       ) {
  K <- x$K
  var.names <- x$var.names
  pars <- x$pars
  dimG <- length(var.names$Gnames)
  dimZ <- length(var.names$Znames)
  valueGtoX <- as.vector(t(x$pars$beta[, -1]))
  valueXtoZ <- as.vector(t(x$pars$mu))
  valueXtoY <- as.vector(x$pars$gamma$beta)[1:K]
  GtoX <- data.frame(source = rep(x$var.names$Gnames, K),
                     target = paste0("Latent Cluster", as.vector(sapply(1:K, function(x) rep(x, dimG)))),
                     value = abs(valueGtoX),
                     group = as.factor(valueGtoX > 0))
  XtoZ <- data.frame(source = paste0("Latent Cluster", as.vector(sapply(1:K, function(x) rep(x, dimZ)))),
                     target = rep(var.names$Znames, K),
                     value = abs(valueXtoZ),
                     group = as.factor(valueXtoZ > 0))
  XtoY <- data.frame(source = paste0("Latent Cluster", 1:K),
                     target = rep(var.names$Ynames, K),
                     value = abs(valueXtoY),
                     group = as.factor(valueXtoY > 0))
  # if(x$family == "binary"){
  #   XtoY$value <- exp(valueXtoY)
  # }
  links <- rbind(GtoX, XtoZ, XtoY)
  nodes <- data.frame(name = unique(c(as.character(links$source), as.character(links$target))),
                      group = as.factor(c(rep("exposure", dimG), 
                                          rep("lc", K), 
                                          rep("biomarker", dimZ), "outcome")))
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  color_scale <- data.frame(domain = c("exposure", "lc", "biomarker", "outcome", "TRUE", "FALSE"),
                            range = c(G_color, X_color, Z_color, Y_color, pos_link_color, neg_link_color))
  
  p <- sankeyNetwork(Links = links, 
                     Nodes = nodes,
                     Source = "IDsource", 
                     Target = "IDtarget",
                     Value = "value", 
                     NodeID = "name",
                     colourScale = JS(
                       sprintf(
                       'd3.scaleOrdinal()
                        .domain(%s)
                        .range(%s)
                       ',
                       jsonlite::toJSON(color_scale$domain),
                       jsonlite::toJSON(color_scale$range)
                     )), 
                     LinkGroup ="group", 
                     NodeGroup ="group",
                     sinksRight = FALSE, 
                     fontSize = fontsize)
  p
}