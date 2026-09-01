## Test environments

* local Fedora 44, R 4.6.0
* GitHub Actions: macOS, Windows, Ubuntu (release and devel)
* R-hub

## R CMD check results

One NOTE, on every platform:

```
* checking foreign function calls ... NOTE
Registration problems:
  Evaluating 'sym$name' during check gives error 'object 'sym' not found':
   .C(sym$name, ..., PACKAGE = sym$dll)
   .Call(sym$name, ..., PACKAGE = sym$dll)
```

This is inherent to what the package does. cppDE generates C++ for the model a
user defines, compiles it at run time and calls the resulting entry point. That
entry point lives in a shared object that does not exist when cppDE is
installed, so there is nothing to register and no symbol a static check can
resolve. The name and the shared object it belongs to are both determined at
run time, which is what `sym$name` and `sym$dll` carry.

The lookup is scoped to the shared object the model was compiled into, so a
model can only reach its own entry points, and a missing symbol is reported by
name rather than dereferenced.

## Downstream dependencies

None on CRAN.
