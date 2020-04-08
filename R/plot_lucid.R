plot.lucid <- function(x){
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
                     colourScale=my_color, LinkGroup="group", NodeGroup="group",
                     sinksRight=FALSE, fontSize = 7)
  p
}