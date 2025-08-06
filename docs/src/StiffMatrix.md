# Stiffness Matrix

_The stiffness matrix_ represents the discretized weak form of differential operators (e.g., curl-curl, grad-div, Laplacian) in finite element methods.

- It encodes the action of the differential operator on the chosen basis functions.
- The stiffness matrix is typically sparse and symmetric (for self-adjoint operators).
- Efficient assembly and application are crucial for large-scale simulations.

## StiffMatrix in FEMSei

The `StiffMatrix` function assembles the global stiffness matrix for various finite element types (e.g., `HCurlElement`, `HDivElement`, `ScalarElement`, `VectorElement`).  
It uses numerical quadrature and the mapping from reference to physical elements.

### General 
```markdown
```@autodocs
CGDycore.FEMSei.VortCrossVel!
```
### Notes

- The function assumes that basis functions and their derivatives are provided in the element objects.
- The Jacobian and its inverse are used for mapping and integration.
- Only nonzero entries are stored in the sparse.

