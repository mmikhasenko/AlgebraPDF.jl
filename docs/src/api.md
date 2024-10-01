# API Documentation

The API for the package

## Standard functions

There is a range of predefined functions
implemented as a struct that holds tuple of parameters

### Parameters

The default parameter type is the `NamedTuple`. It holds parameter values and names and can update their value.
The additional type `FlaggedNamedTuple` allows to group parameters to `freed` and `fixed`

```@docs
FlaggedNamedTuple
```

The package introduces two operations on parameter types, `merge` and `subtract`.
The `merge(nt1::NamedTuple, nt1::NamedTuple)` is defined in `Base`.

```@docs
subtract
```

### Functions

Several common function are defined

```@docs
FGauss
FExp
FPol
FDoubleGaussFixedRatio
```

```@docs
AlgebraPDF.noparsnormf
AlgebraPDF.nt
```

### Macros

```@docs
AlgebraPDF.@makepdftype
AlgebraPDF.@makefuntype
```

### Other Functions

```@docs
AlgebraPDF.updatevalueorflag
```
