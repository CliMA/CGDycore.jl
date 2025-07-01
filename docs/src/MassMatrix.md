# Massmatrix

_The mass matrix_ represents the inner product in the discretized function space.

 - For efficiency, we avoid explicitly assembling the mass matrix and instead
  apply its action directly to fields.
 - In Julia, this is often achieved using **matrix-free** methods and
  broadcasting, enabling efficient and flexible computations.

The mass matrix acts as a "pseudo-operator": it cannot be called directly, but
can be applied to fields in a manner similar to other operators, leveraging
broadcasting for performance and composability.

## MassMatrix with HCurlElement 
```@docs
CGDycore.FEMSei.MassMatrix
```