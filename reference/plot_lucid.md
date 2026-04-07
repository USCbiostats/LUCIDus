# Visualize LUCID model through a Sankey diagram

In the Sankey diagram, each node either represents a variable (exposure,
omics or outcome) or a latent cluster. Each line represents an
association. The color of the node represents variable type, either
exposure, omics or outcome. The width of the line represents the effect
size of a certain association; the color of the line represents the
direction of a certain association.

## Usage

``` r
plot_lucid(
  x,
  G_color = "dimgray",
  X_color = "#eb8c30",
  Z_color = "#2fa4da",
  Y_color = "#afa58e",
  pos_link_color = "#67928b",
  neg_link_color = "#d1e5eb",
  fontsize = 7
)
```

## Arguments

- x:

  A LUCID model fitted by [`est_lucid`](est_lucid.md)

- G_color:

  Color of node for exposure

- X_color:

  Color of node for latent cluster

- Z_color:

  Color of node for omics data

- Y_color:

  Color of node for outcome

- pos_link_color:

  Color of link corresponds to positive association

- neg_link_color:

  Color of link corresponds to negative association

- fontsize:

  Font size for annotation

## Value

A DAG graph created by `sankeyNetwork`

## Examples

``` r
if (FALSE) { # \dontrun{
# prepare data
G <- sim_data$G
Z <- sim_data$Z
Y_normal <- sim_data$Y_normal
Y_binary <- sim_data$Y_binary
cov <- sim_data$Covariate

# plot lucid model
fit1 <- est_lucid(G = G, Z = Z, Y = Y_normal, CoY = NULL, family = "normal", 
K = 2, seed = 1008)
plot_lucid(fit1)

# change node color
plot_lucid(fit1, G_color = "yellow")
plot_lucid(fit1, Z_color = "red")

# change link color
plot_lucid(fit1, pos_link_color = "red", neg_link_color = "green")
} # }
```
