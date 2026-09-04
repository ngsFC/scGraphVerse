# Toy adjacency matrix for examples

An adjacency matrix generated using real dataset and STRINGdb for
demonstrating functions in the scGraphVerse package.

## Usage

``` r
data(toy_adj_matrix)
```

## Format

A matrix of dimension 10x10.

## Examples

``` r
data(toy_adj_matrix)
str(toy_adj_matrix)
#>  num [1:35, 1:35] 0 1 1 0 0 0 0 1 0 0 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:35] "ACTG1" "ARPC2" "ARPC3" "BTF3" ...
#>   ..$ : chr [1:35] "ACTG1" "ARPC2" "ARPC3" "BTF3" ...
```
