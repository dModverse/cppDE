# Reject Symbol Names the Backend Cannot Parse

Stops when an expression or a symbol name holds a Python keyword.

## Usage

``` r
checkSymbolNames(...)
```

## Arguments

- ...:

  Character vectors of expressions or symbol names, optionally named.
  Names are checked along with the values. `NULL` is ignored.

## Value

`NULL`, invisibly. Called for the error.
