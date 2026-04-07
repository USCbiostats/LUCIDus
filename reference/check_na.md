# Check missing patterns in omics data Z

Check missing patterns in omics data Z

## Usage

``` r
check_na(Z)
```

## Arguments

- Z:

  A data matrix representing omics data

## Value

1\. index:indeces for missing values in omics data 2. indicator_na:
missing pattern for each observation 3. impute_flag: - flag to
initialize imputation. Only happens when sporadic missing pattern is
observed
