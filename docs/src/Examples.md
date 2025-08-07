# Examples

This section provides various example initializations for numerical experiments with CGDycore.

## Bickley Jet Planar

The classic Bickley jet setup is a standard test case for barotropic instability and nonlinear jet evolution in atmospheric and oceanic flows.
```@docs
CGDycore.Examples.BickleyJetExample
```

## Modon Collision

The Modon collision example initializes two interacting modons (coherent vortex structures), which are used to study vortex dynamics, wave interactions, and nonlinear evolution in geophysical flows.
```@docs
CGDycore.Examples.ModonCollisionExample
```

## Galewsky

The Galewsky test case is a standard barotropic instability benchmark on the sphere, featuring a zonal jet and a localized perturbation. It is widely used to assess the performance of global shallow water models.
```@docs
CGDycore.Examples.GalewskyExample
```

## Rossby-Haurwitz Wave

The Haurwitz wave example initializes a steady, analytical wave solution on the sphere. It is commonly used to test the advection and wave propagation properties of numerical models.
```@docs
CGDycore.Examples.HaurwitzExample
```

## More Examples

Additional initializations and test cases are included in the `Examples` module and can be used analogously.
