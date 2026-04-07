# generate bootstrp ci (normal, basic and percentile)

generate bootstrp ci (normal, basic and percentile)

## Usage

``` r
gen_ci(x, conf = 0.95)
```

## Arguments

- x:

  an object return by boot function

- conf:

  A numeric scalar between 0 and 1 to specify confidence level(s) of the
  required interval(s).

## Value

a matrix, the first column is t0 statistic from original model
