#' Visualize the LUCID model through a Sankey diagram
#' This function generates a Sankey diagram for the results of integrative clustering based on an \code{lucid} object
#' @param x A model fitted by \code{\link{est.lucid}}
#' @param ... Other parameters to be passed to \code{plot}
#' @return A DAG graph created by \code{\link{sankeyNetwork}}
#' @import networkD3
#' @export
#' @author Cheng Peng, Zhao Yang, David V. Conti
#' @references
#' Cheng Peng, Jun Wang, Isaac Asante, Stan Louie, Ran Jin, Lida Chatzi, Graham Casey, Duncan C Thomas, David V Conti, A Latent Unknown Clustering Integrating Multi-Omics Data (LUCID) with Phenotypic Traits, Bioinformatics, , btz667, https://doi.org/10.1093/bioinformatics/btz667.
#'
#' @examples
#' \dontrun{
#' fit1 <- est.lucid(G = G1, Z = Z1, Y = Y1, CoY = CovY, K = 2, family = "binary")
#' plot(fit1)
#' }
plot.lucid <- function(x, ...){
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
  if(x$family == "binary"){
    XtoY$value <- exp(valueXtoY)
  }
  links <- rbind(GtoX, XtoZ, XtoY)
  nodes <- data.frame(name = unique(c(as.character(links$source), as.character(links$target))),
                      group = as.factor(c(rep("exposure", dimG), 
                                          rep("lc", K), 
                                          rep("biomarker", dimZ), "outcome")))
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  my_color <- 'd3.scaleOrdinal() .domain(["exposure", "lc", "biomarker", "outcome", "TRUE", "FALSE"]) .range(["dimgray", "#eb8c30", "#2fa4da", "#afa58e", "#67928b", "#d1e5eb"])'
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name",
                     colourScale = my_color, LinkGroup ="group", NodeGroup ="group",
                     sinksRight = FALSE, fontSize = 7)
  p
}