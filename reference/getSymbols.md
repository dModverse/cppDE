# Extract Symbol Names from R Expressions

Returns the unique names of all symbols occurring in a character vector
of R expressions. `pi` is the number and no symbol.

## Usage

``` r
getSymbols(expr, omit = NULL)
```

## Arguments

- expr:

  Character vector of R expressions.

- omit:

  Optional character vector of symbol names to remove.

## Value

Character vector of unique symbol names.
