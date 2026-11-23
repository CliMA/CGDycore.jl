# The continuous Galerkin spectral element dynamical core of CGDycore.jl (`src/CGSEM`)

*Numerical scheme and comparison with HOMME and ClimaCore*

Summary of the implementation

*Converted from `CGSEMScheme.tex`; the PDF is the reference for layout.*

**Abstract.** This note describes the numerical scheme implemented in the module
`CGSEM` of CGDycore.jl as it is coded in the kernels of
`src/CGSEM`. The scheme combines a continuous Galerkin spectral element
method (SEM) on Gauss–Lobatto–Legendre (GLL) points in the horizontal with
a low-order, staggered (Charney–Phillips type) discretization in a
terrain-following height coordinate in the vertical. The prognostic variables
are density, the horizontal velocity components and the vertical velocity
(on layer interfaces) in the local orthonormal spherical frame
(east, north, radial), the density-weighted potential temperature and tracer
densities. The momentum equation is used in vector-invariant form, and all
metric terms of the sphere enter through the adjugate of the element mapping
written in the spherical frame. Time integration is done by Rosenbrock-type W-methods with an
approximate Jacobian that contains only the vertical acoustic and buoyancy
coupling (HEVI). We describe the discrete operators term by term and compare
the design with the spectral element cores HOMME (E3SM) and ClimaCore
(CliMA).

## 1 Overview

The main routine for the right-hand side is `Fcn!` in
`CGSEM/Fcn.jl`. One evaluation consists of

1. a diagnostic step `FcnPrepare!` (pressure, temperature, eddy
  coefficients),
2. a first pass of horizontal hyperviscosity (Laplacians of momentum,
  potential temperature and tracers), followed by direct stiffness summation
  (DSS) and parallel exchange,
3. the second pass, in which all tendencies are computed: biharmonic
  hyperviscosity, rotational and Coriolis terms, gradient of kinetic energy,
  pressure gradient and gravity, transport of $\rho$, $\rho\theta$ and the
  tracers,
4. DSS/parallel exchange of the tendency and the column physics (vertical
  diffusion, surface fluxes, forcing, microphysics, sedimentation, damping).

Before the first pass the meridional velocity is set to zero at the
boundary degrees of freedom `BoundaryDoF` (lateral walls of Cartesian
channel configurations; the set is empty on the sphere). If the tracer
limiter is active, the element bounds of the tracers are computed and their
exchange is started at the same time.

All element-wise operators are evaluated in KernelAbstractions kernels, one
thread per horizontal node and layer; a work group holds the $N\times N$
nodes of one element and a tile of layers, and the neighbouring layers that
a vertical stencil needs are loaded into shared memory as halo values. The contributions of the elements are
accumulated into the global degrees of freedom by atomic additions and divided
by the assembled lumped mass. In this way DSS is part of every operator
kernel. The elements are ordered so that the elements at the boundary of a
process come first: they are computed first, the exchange is started, the
interior elements are computed and finally the exchange is completed (overlap
of communication and computation).

**Table 1:** Main components of the CGSEM scheme.

| Horizontal | CG spectral elements on quadrilaterals (cubed sphere or Cartesian), GLL collocation, polynomial order $p$ (typically $p=3$ or $4$), lumped mass, DSS |
|---|---|
| Vertical | terrain-following height coordinate, $N_z$ layers, lowest-order staggering: $\rho,u,v,\rho\theta,\rho q$ in layers, $w$ on interfaces |
| Momentum | vector-invariant form, prognostic velocities $u,v,w$ in the local spherical frame $(\boldsymbol e_\lambda,\boldsymbol e_\varphi,\boldsymbol e_r)$ |
| Metric | adjugate of the element mapping in the spherical frame, stored at the $N\times N\times 2$ nodes of every cell |
| Thermodynamics | $\rho\theta$ (dry); variants with internal or total energy |
| Horizontal transport | centred SEM operators; split form for $\rho\theta$ |
| Vertical transport | limited third-order upwind reconstruction (Koren-type limiter, non-uniform grid) |
| Dissipation | biharmonic hyperviscosity, separate coefficients for vorticity, divergence and scalars |
| Tracer limiter | optional optimization-based limiter (local bounds, mass conservation) |
| Time integration | Rosenbrock W-methods (e.g. `SSP-Knoth`) with a frozen vertical acoustic Jacobian, Schur complement for $w$, tridiagonal solves per column |

## 2 Continuous equations

The dry compressible Euler equations are solved in the vector-invariant form

$$
\begin{aligned}
  \partial_t \rho + \nabla\cdot(\rho\boldsymbol{v}) &= 0, \qquad (1) \\
  \partial_t \boldsymbol{v} + (\boldsymbol{\omega} + 2\boldsymbol{\Omega})\times\boldsymbol{v}
    + \nabla K + \frac{1}{\rho}\nabla p + \nabla\Phi &= \boldsymbol{D}_{\boldsymbol{v}}, \qquad (2) \\
  \partial_t (\rho\theta) + \nabla\cdot(\rho\theta\boldsymbol{v}) &= D_{\theta}, \qquad (3) \\
  \partial_t (\rho q_l) + \nabla\cdot(\rho q_l\boldsymbol{v}) &= D_{q_l} + S_{q_l},
    \qquad l = 1,\dots,N_{Tr}, \qquad (4)
\end{aligned}
$$

with $\boldsymbol{\omega} = \nabla\times\boldsymbol{v}$, the kinetic energy
$K = \tfrac12|\boldsymbol{v}|^2$, the geopotential $\Phi$ (gravity
$\boldsymbol g = -\nabla\Phi$, shallow or deep), the dissipation terms $\boldsymbol D$ and
the physical sources $S$. For the dry state the pressure is

$$
p = p_0\left(\frac{R_d\,\rho\theta}{p_0}\right)^{1/(1-\kappa)},
  \qquad \kappa = R_d/c_{pd}. \tag{5}
$$

On the sphere the velocity is used in its spherical representation,

$$
\boldsymbol{v} = u\,\boldsymbol e_\lambda + v\,\boldsymbol e_\varphi + w\,\boldsymbol e_r , \tag{6}
$$

i.e. $u$, $v$ and $w$ are the components with respect to the local
orthonormal basis of the spherical coordinates $(\lambda,\varphi,r)$ (east,
north, radial) at the node, not Cartesian components of $\mathbb R^3$. The
basis rotates from node to node. The formulation is derived from the
Cartesian vector-invariant form by moving the rotation matrix into the
contravariant metric (Section 3.3), so that the rotation is
not represented explicitly in the operators. In the Cartesian configurations (channel,
box) $u$, $v$, $w$ are the Cartesian components.

The Coriolis term and gravity are prescribed by the function objects
`CoriolisFun` and `GravitationFun`. The variants that fit the
CGSEM kernel interface are `CoriolisShallow`
($2\boldsymbol\Omega\times\boldsymbol{v}$ reduced to the vertical component
$f = 2\Omega\sin\varphi$), `FPlane` ($f = f_0 + \beta(y-y_0)$),
`CoriolisNo`, which act only on $u$ and $v$, and
`CoriolisNonTraditional` (`--CoriolisType="NonTraditional"` or
`"Deep"`), the full $2\boldsymbol\Omega\times\boldsymbol{v}$ including the
non-traditional terms with $f_c = 2\Omega\cos\varphi$. Gravity is
$g$ (`GravitationShallow`) or $g\,(a/r)^2$ (`GravitationDeep`),
evaluated at the node.

Other options select internal or total energy instead of
$\rho\theta$ (`DryInternalEnergy`, `DryTotalEnergy`, …), moisture
variables, a turbulent kinetic energy and EDMF variables; they reuse the same
operators and are not discussed here.

## 3 Grid, coordinate transformation and metric terms

### 3.1 Horizontal elements and vertical layers

The horizontal mesh consists of quadrilateral elements $\Omega_e$, for
example the elements of an equiangular or equidistant cubed sphere with
$n_{\text{Panel}}\times n_{\text{Panel}}$ elements per panel. Every column of
elements is divided into $N_z$ layers between the interfaces
$z_{1/2} < z_{3/2} < \dots < z_{N_z+1/2} = H$ (optionally stretched). A
grid cell $\Omega_{e,k}$ is mapped from the reference cell
$(\xi^1,\xi^2,\zeta)\in[-1,1]^3$ by

$$
\boldsymbol{x}(\xi^1,\xi^2,\zeta) = \frac{\boldsymbol{X}(\xi^1,\xi^2)}{|\boldsymbol{X}(\xi^1,\xi^2)|}
    \bigl(a + h(z(\zeta), z_s(\xi^1,\xi^2))\bigr),
  \quad
  z(\zeta) = \tfrac12\bigl((1-\zeta) z_{k-1/2} + (1+\zeta) z_{k+1/2}\bigr), \tag{7}
$$

where $\boldsymbol X$ is the bilinear interpolation of the four element corners on the
cube, $a$ is the radius of the planet and $h(z,z_s)$ is the terrain-following
transformation (`AdaptGrid`, e.g. Gal-Chen or SLEVE) of the
orography $z_s$. In the Cartesian case $\boldsymbol x$ is the analogous mapping
without the projection onto the sphere.

### 3.2 Nodes and derivative matrices

In the horizontal the GLL points $\xi_i$ with the weights $\omega_i$,
$i=1,\dots,N$, $N=p+1$, are used as interpolation and quadrature nodes. The
strong and weak derivative matrices are

$$
D_{ij} = \ell_j'(\xi_i), \qquad
  D^W = -W^{-1} D^T W, \quad W = \operatorname{diag}(\omega_1,\dots,\omega_N), \tag{8}
$$

(`DS` and `DW` in the code). In the vertical each layer has the
two GLL points $\zeta = -1$ (bottom, index $K=1$) and $\zeta=+1$ (top,
$K=2$); the corresponding two-point quadrature has the weights $1,1$.
The metric terms are stored at the $N\times N\times 2$ points of every cell.

### 3.3 Cartesian contravariant metric and rotation to the spherical frame

The idea of the formulation is that all information about the sphere is
contained in the metric terms. The starting point is the vector-invariant
form (1)–(4) with the Cartesian velocity
$\boldsymbol{v}_C\in\mathbb R^3$ and the Cartesian contravariant metric

$$
J = \det\frac{\partial\boldsymbol x}{\partial\boldsymbol\xi},\qquad
  \mathcal{J}^a_C = J\,\nabla\xi^a \;\in\mathbb R^3,\qquad a = 1,2,3, \tag{9}
$$

($\xi^3=\zeta$), for which the operators

$$
J\nabla\cdot\boldsymbol f = \sum_{a}\partial_{\xi^a}\bigl(\mathcal{J}^a_C\cdot\boldsymbol f\bigr),\quad
  J\nabla\phi = \sum_a \mathcal{J}^a_C\,\partial_{\xi^a}\phi,\quad
  J\nabla\times\boldsymbol{v}_C = \sum_a \partial_{\xi^a}\bigl(\mathcal{J}^a_C\times\boldsymbol{v}_C\bigr) \tag{10}
$$

are exact because of the metric identities $\sum_a\partial_{\xi^a}\mathcal{J}^a_C = 0$.

Let $R(\lambda,\varphi)$ be the rotation matrix whose rows are the unit
vectors $\boldsymbol e_\lambda,\boldsymbol e_\varphi,\boldsymbol e_r$ at the node
(`Grids.MCart2Sphere`). The spherical velocity is

$$
\boldsymbol{v} = R\,\boldsymbol{v}_C,\qquad \boldsymbol{v}_C = R^T\boldsymbol{v},\qquad \boldsymbol{v} = (u,v,w)^T . \tag{11}
$$

Substituting $\boldsymbol{v}_C = R^T\boldsymbol{v}$ into the Cartesian forms, $R$ can be moved
into the metric. The spherical contravariant metric is

$$
\mathcal{J}^a = R\,\mathcal{J}^a_C, \qquad
  \mathcal{J}^a_{\;b} = \sum_c R_{bc}\,(\mathcal{J}^a_C)_c = \mathcal{J}^a_C\cdot\boldsymbol e_b, \tag{12}
$$

i.e. the adjugate of the Cartesian Jacobian multiplied by the rotation
matrix. This is what the code stores in `dXdxI[a,b,K,ID,k,e]`
($=\mathcal{J}^a_{\;b}$ at vertical node $K$) together with `J[ID,K,k,e]`.

**Computation of the metric.**
All metric terms are computed in `FiniteElements/Metric.jl`
(`MetricCompute`).

- The node positions $\boldsymbol X$ are generated
  (`FillXKernel!`), the Cartesian adjugate is computed from the
  derivatives of $\boldsymbol X$ by cross products (single differentiation,
  `--MetricType="Single"`, `MetricQuadKernel!`, default for CG) or
  in the curl form of Kopriva (double differentiation,
  `--MetricType="Double"`, `FillContraCurlQuadKernel1-3!`, default
  for DG), which fulfils the discrete metric identities. With only two
  vertical nodes per layer the double form evaluates $\mathcal{J}^{1,2}$ with the
  mean radius $\bar r$ of the layer instead of the node radius $r_K$, an
  error $\Delta r/(2r)$ that also enters $J$. `ApplyRotateKernel!` then multiplies with
  $R^T$ from the right as in (12), and $J = |\det\mathcal{J}|^{1/2}$
  (`DetKernel!`, for CGSEM `DetCGKernel!` with $J$ stored per
  element node). This is the exact rotated Cartesian (deep) metric of the
  polynomial interpolant of $\boldsymbol x$. For the shallow atmosphere (every
  `Equation` except `CompressibleDeep`) the tangential rows of $L$
  are scaled by $a/r$, $L\to SL$ with $S = \operatorname{diag}(a/r,a/r,1)$,
  i.e. $\mathcal{J}^a_{\;b}\to\mathcal{J}^a_{\;b}\,s_b$ with $s = (a/r,\,a/r,\,(a/r)^2)$ and
  $J\to(a/r)^2J$ (`ShallowMetricKernel!`). This is the shallow metric of
  Guba below, evaluated with the polynomial derivatives of $\boldsymbol x$ instead of
  the analytic ones; without orography the two agree to $2\cdot10^{-7}$
  ($n_{\text{Panel}} = 16$, $p = 3$) and converge with order $p+2$.
- The earlier metric of Guba (`JacobiSphere3Kernel!`, removed from
  `Grids`) has the same continuous form. Here the Jacobian is formed directly in the spherical frame as a product
  of two matrices: the derivative $\tfrac14 BC = \partial\boldsymbol X/\partial\boldsymbol\xi$
  of the bilinear cube point and the matrix
  $\tfrac{a}{|\boldsymbol X|}\,D_\lambda A$ (`DD*A`) that projects onto the sphere and onto
  $\boldsymbol e_\lambda$, $\boldsymbol e_\varphi$,

  $$
L = \begin{pmatrix}
      \boldsymbol e_\lambda\cdot\partial_{\xi^1}\boldsymbol x_s & \boldsymbol e_\lambda\cdot\partial_{\xi^2}\boldsymbol x_s & 0\\
      \boldsymbol e_\varphi\cdot\partial_{\xi^1}\boldsymbol x_s & \boldsymbol e_\varphi\cdot\partial_{\xi^2}\boldsymbol x_s & 0\\
      \partial_{\xi^1} h & \partial_{\xi^2} h & \partial_\zeta h
    \end{pmatrix},\qquad \boldsymbol x_s = a\frac{\boldsymbol X}{|\boldsymbol X|}, \tag{13}
$$

  and $\mathcal{J} = \operatorname{adj}L$, $J = \det L$. Since the rows of $R$ are
  orthonormal, $L = R\,\partial\boldsymbol x/\partial\boldsymbol\xi$ and
  $\operatorname{adj}L = \operatorname{adj}(\partial\boldsymbol x/\partial\boldsymbol\xi)\,R^T$,
  which is (12). The shallow-atmosphere version uses the radius
  $a$ in the tangential rows instead of $a+h$; the deep version scales them
  by $(a+h)/a$ and is the exact rotated Cartesian metric. The entries $L_{13} = L_{23} = 0$ reflect that
  the $\zeta$-lines are radial. The derivatives of $h$ are computed with $D$
  in the horizontal and by
  $\partial_\zeta h \approx \tfrac12(h_{K=2}-h_{K=1})$ in the vertical.

In the Cartesian configurations $R = I$ (`ApplyRotateKernel!` is
empty).

**Metric in reduced precision.**
The metric kernels differentiate the element-local positions
$\boldsymbol X - \boldsymbol X_1$ ($\boldsymbol X_1$ the first node of the element). These are formed
without the absolute positions as
$(h-h_1)\,\hat{\boldsymbol e} + (a+h_1)(\hat{\boldsymbol e}-\hat{\boldsymbol e}_1)$, with $\hat{\boldsymbol e}$ the
unit radial direction and $h$ the height of the node, so that vertical
differences are exact height differences (`XLocal!`). Below Float64 the
positions, the contravariant metric and $\Delta z$ are computed in Float64 on the
CPU and then rounded to the working precision (`OnMetricBackend!`); this
also works on backends without Float64 such as Metal. In Float32,
$|\boldsymbol X|\approx 6.4\cdot10^6\,$m is rounded to about $0.5\,$m, which is a large
fraction of the lowest layers, and in the curl form the products
$X_l\,\partial X_m$ (element size squared) cancel to $\mathcal{J}^{1,2}$ (element size
times layer depth). For CGSEM the rotation (or the normal column and the tangential split of the
Cartesian frame), $J$ and the shallow scaling are part of this Float64 setup
(`OnFrameBackend!`): rotated in Float32, the large vertical row $\mathcal{J}^3$ left
a remainder $\epsilon|\mathcal{J}^3|$ in its tangential components, a spurious horizontal
force $\epsilon g$. Before, $J$ had relative errors of $5\cdot10^{-4}$ (CG) and
$7\cdot10^{-3}$ (DG, curl form), the DG face terms up to $2\cdot10^{-2}$ and the
tangential components of $\mathcal{J}^3$ $6\cdot10^{-7}|\mathcal{J}^3|$; now all metric terms
agree with Float64 to $10^{-6}$ and the tangential components of $\mathcal{J}^3$ to
$7\cdot10^{-9}|\mathcal{J}^3|$. In the baroclinic wave ($n_{\text{Panel}}=32$, 64 levels,
two days) the Float32 solution of the spherical frame then deviates from the
Float64 one by $1.2\cdot10^{-5}$ in $\rho$, $3.4\cdot10^{-3}$ in $v$ and $2.4\cdot10^{-2}$
in $w$ (before $3\cdot10^{-4}$, $9\cdot10^{-2}$ and $0.95$).

**Lam\'e form.**
Equivalently, $L = \operatorname{diag}(h_\lambda,h_\varphi,h_r)\, \partial(\lambda,\varphi,r)/\partial\boldsymbol\xi$ with $h_\lambda = r\cos\varphi$,
$h_\varphi = r$, $h_r = 1$ (shallow: $r\to a$). With $H = h_\lambda h_\varphi h_r$
and $\hat J = \det\partial(\lambda,\varphi,r)/\partial\boldsymbol\xi$,

$$
J = H\hat J,\qquad
  \mathcal{J}^a_{\;b} = \frac{H}{h_b}\,\hat J\,\frac{\partial\xi^a}{\partial q^b},
  \qquad
  \sum_a\partial_{\xi^a}\bigl(\mathcal{J}^a_{\;b}\,f\bigr)
   = \hat J\,\partial_{q^b}\Bigl(\frac{H}{h_b}\,f\Bigr), \tag{14}
$$

$(q^1,q^2,q^3) = (\lambda,\varphi,r)$. The spherical metric does not fulfil
$\sum_a\partial_{\xi^a}\mathcal{J}^a_{\;b} = 0$; the non-zero remainder is
$\hat J\,\partial_{q^b}(H/h_b)$, i.e. exactly the curvature of the basis
vectors.

### 3.4 Discrete differential operators in spherical components

The code uses the operators (10) with $\mathcal{J}^a_C$, $\boldsymbol{v}_C$
replaced by $\mathcal{J}^a$, $\boldsymbol{v}$:

$$
J\nabla\cdot\boldsymbol f = \sum_{a}\partial_{\xi^a}\bigl(\mathcal{J}^a\cdot\boldsymbol f\bigr),\quad
  J\nabla\phi = \sum_a \mathcal{J}^a\,\partial_{\xi^a}\phi,\quad
  J\boldsymbol\omega = \sum_a \partial_{\xi^a}\bigl(\mathcal{J}^a\times\boldsymbol{v}\bigr). \tag{15}
$$

Whether this is the Cartesian scheme written in spherical components depends
on whether the rotation $R$ appears inside or outside a derivative.

- *Fluxes and divergence.*
  $\mathcal{J}^a\cdot\boldsymbol{v} = (R\mathcal{J}^a_C)\cdot(R\boldsymbol{v}_C) = \mathcal{J}^a_C\cdot\boldsymbol{v}_C$ holds
  pointwise. The contravariant mass, $\rho\theta$ and tracer fluxes, and
  hence the discrete continuity, $\rho\theta$ and tracer equations, are
  identical node by node to the Cartesian discretization. The same applies to
  the lower boundary condition $\mathcal{J}^3\cdot\boldsymbol{v} = 0$.
- *Gradient.*
  $\sum_a\mathcal{J}^a\partial_{\xi^a}\phi = R\sum_a\mathcal{J}^a_C\partial_{\xi^a}\phi$:
  the discrete gradient (pressure gradient, gradient of $K$) is the rotated
  Cartesian discrete gradient, exactly.
- *Kinetic energy.* $|\boldsymbol{v}| = |\boldsymbol{v}_C|$, so $K$ is identical.
- *DSS.* $R$ depends only on the position of the node, so all
  copies of a node have the same $R$ and the assembly commutes with the
  rotation.
- *Curl.* Here the rotation is inside the derivative. With
  $R(\boldsymbol a\times\boldsymbol b) = (R\boldsymbol a)\times(R\boldsymbol b)$,

  $$
R\sum_a\partial_{\xi^a}\bigl(\mathcal{J}^a_C\times\boldsymbol{v}_C\bigr)
    = \sum_a\partial_{\xi^a}\bigl(\mathcal{J}^a\times\boldsymbol{v}\bigr)
     - \sum_a\bigl(\partial_{\xi^a}R\bigr)R^T\bigl(\mathcal{J}^a\times\boldsymbol{v}\bigr). \tag{16}
$$

  The skew-symmetric matrices $(\partial_{\xi^a}R)R^T$ describe the rotation
  of the local basis from node to node. They are not contained in the
  metric (12). Since the $\zeta$-lines are radial,
  $\partial_\zeta R = 0$ and only the horizontal derivatives contribute. If
  only the first term on the right-hand side is used, the continuous limit
  of (14) gives

  $$
\boldsymbol\omega^{\text{code}} = \nabla\times\boldsymbol{v} + \delta\boldsymbol\omega,\qquad
    \delta\omega_i = \sum_{b,c}\varepsilon_{ibc}\,
      \frac{v_c}{h_i h_b}\,\partial_{q^b}h_i , \tag{17}
$$

  which is

  $$
\delta\boldsymbol\omega = \Bigl(-\frac{v + w\tan\varphi}{r},\;\frac{u}{r},\;0\Bigr)
    \ \ \text{(deep, rotated Cartesian metric)},\qquad
    \delta\boldsymbol\omega = \Bigl(-\frac{w\tan\varphi}{a},\;0,\;0\Bigr)
    \ \ \text{(shallow metric of Guba)}. \tag{18}
$$

  (The deep result was verified symbolically from (16).)
  The vertical vorticity, including the curvature term $u\tan\varphi/r$, is
  reproduced exactly in both cases, because $h_r = 1$. The horizontal
  components carry the additional terms; in the shallow case the momentum
  equation would get $-w^2\tan\varphi/a$ in $\partial_t v$ and
  $+vw\tan\varphi/a$ in $\partial_t w$. These terms are singular at the
  poles, where the cubed sphere with even $n_{\text{Panel}}$ has nodes, so
  they cannot be compensated by analytic source terms.

**Curl of the Cartesian velocity.**
The code therefore evaluates the horizontal derivatives in the curl with the
Cartesian components inside the derivative and rotates the result back,

$$
J\boldsymbol\omega = R\sum_{a=1,2} D_a\bigl(R^T(\mathcal{J}^a\times\boldsymbol{v})\bigr)
     + D_\zeta\bigl(\mathcal{J}^3\times\boldsymbol{v}\bigr), \tag{19}
$$

where $R$ is evaluated at every node from its position (`CurlFrame`),
$R^T(\mathcal{J}^a\times\boldsymbol{v}) = \mathcal{J}^a_C\times\boldsymbol{v}_C$, and the vertical part is
unchanged because $R$ is constant along a column. With the deep metric this
is node by node the rotated Cartesian curl, also at the poles. With the
shallow metric, (14) applied to (19) gives in the
continuous limit $\nabla\times\boldsymbol{v} + (v/a,\,-u/a,\,0)$, i.e. the additional
tendencies $(uw/a,\; vw/a,\; -(u^2+v^2)/a)$, which contain no
$\tan\varphi$. They are subtracted in every half cell,

$$
\partial_t u \mathrel{-}= \frac{u w}{a},\qquad
  \partial_t v \mathrel{-}= \frac{v w}{a},\qquad
  \partial_t w \mathrel{+}= \frac{u^2 + v^2}{a}, \tag{20}
$$

with the half-cell weights $J^{(K)}$ and the half-cell vertical velocities
$w_{k\mp1/2}$ of the rotational term (`CurvatureCoefficient` $=1/a$, zero
for `CompressibleDeep` and Cartesian grids). The correction is energy
neutral in each half cell, $u\,(uw) + v\,(vw) - w\,(u^2+v^2) = 0$. For
$u = U\cos\varphi$, $v = 0$, $w = W$ the discrete rotational term converges
with order $p$ to $-\boldsymbol\omega\times\boldsymbol{v}$ of the shallow and of the deep
equations, with Guba's and with the `DGMetric` metric, also on grids
with nodes at the poles.

### 3.5 Lumped mass matrix and direct stiffness summation

Let $g$ be a global horizontal node and $k$ a layer. The assembled half-cell
masses are

$$
M^{(K)}_{g,k} = \sum_{(e,i,j)\mapsto g} J^{(K)}_{e,ij,k}, \qquad K=1,2,
  \qquad M_{g,k} = M^{(1)}_{g,k} + M^{(2)}_{g,k} \tag{21}
$$

(`MassCGKernel!`, array `M[k,g,K]`). The sum runs over all copies of
the node $g$ in the elements that contain it. The exact lumped GLL mass would
be $\omega_i\omega_j M_{g,k}$; on a conforming quadrilateral mesh the product
$\omega_i\omega_j$ is the same for all copies of a node (the GLL weights are
symmetric), so the weights cancel in every quotient
``assembled contribution / assembled mass” and are omitted in both. The
weighted mass $\sum\tfrac12(J^{(1)}+J^{(2)})\,\omega_i\omega_j$ is assembled
separately (`MMass`) and is used only where real volumes are needed,
e.g. for global integrals.

A generic term is evaluated in strong form element by element and assembled
as a $J$-weighted average,

$$
F_{g,k} = \frac{1}{M_{g,k}} \sum_{(e,i,j)\mapsto g} R_{e,ij,k}, \tag{22}
$$

where $R_{e,ij,k}$ is $J$ times the pointwise tendency, summed over the two
vertical nodes of the layer. This is the SEM collocation with lumped mass
together with DSS.
The vertical velocity lives on the interfaces $k+\tfrac12$ and gets the
assembled interface masses $M^{(2)}_{g,k}+M^{(1)}_{g,k+1}$, weighted with the
density (Section 5.6).

## 4 Vertical staggering and boundary conditions

Per horizontal node and layer $k=1,\dots,N_z$ the unknowns
$\rho_k, u_k, v_k, (\rho\theta)_k, (\rho q_l)_k$ are constant in the layer.
The vertical velocity $w_{k+1/2}$, $k=1,\dots,N_z-1$, is stored on the interior
interfaces (array position $k$ belongs to the top of layer $k$). Inside a
layer the velocity at the two vertical nodes is

$$
\boldsymbol{v}^{(1)}_k = (u_k, v_k, w_{k-1/2}), \qquad \boldsymbol{v}^{(2)}_k = (u_k, v_k, w_{k+1/2}), \tag{23}
$$

so $u,v$ are piecewise constant and $w$ is piecewise linear in $\zeta$
(Charney–Phillips staggering of $w$ against $\rho,\rho\theta$).
At the top $w_{N_z+1/2}=0$. At the lower boundary the kinematic condition
$\mathcal{J}^3\cdot\boldsymbol{v} = 0$ (no flow through the orography) gives

$$
w_{1/2} = -\frac{\mathcal{J}^3_{\;1}u_1 + \mathcal{J}^3_{\;2}v_1}{\mathcal{J}^3_{\;3}}. \tag{24}
$$

The vertical fluxes of mass, $\rho\theta$ and tracers vanish at the bottom and
at the top.

## 5 Spatial discretization of the individual terms

In the following, sums over $m$ are horizontal SEM derivatives inside an
element, e.g. $(D_1 f)_{ij} = \sum_m D_{im}f_{mj}$ and
$(D_2 f)_{ij} = \sum_m D_{jm}f_{im}$. Every element contribution is assembled
with (22).

### 5.1 Continuity equation

The horizontal mass fluxes in contravariant form, summed over both vertical
nodes of the layer, are

$$
U^a_k = \rho_k \sum_{K=1}^2 \bigl(\mathcal{J}^{a,(K)}_{\;1} u_k + \mathcal{J}^{a,(K)}_{\;2} v_k\bigr),
  \qquad a = 1,2, \tag{25}
$$

(`Contra12`). The vertical mass flux through the interface $k+\tfrac12$
uses the metric from both sides,

$$
W_{k+1/2} = \rho_k\,\mathcal{J}^{3,(2)}_k\cdot(u_k,v_k,w_{k+1/2})
    + \rho_{k+1}\,\mathcal{J}^{3,(1)}_{k+1}\cdot(u_{k+1},v_{k+1},w_{k+1/2}) \tag{26}
$$

(`Contra3`). The element contribution to the density tendency is

$$
R^\rho_k = -\bigl(D_1 U^1_k + D_2 U^2_k\bigr)
     - \tfrac12\bigl(W_{k+1/2} - W_{k-1/2}\bigr). \tag{27}
$$

### 5.2 Potential temperature: horizontal split form

With $\theta_k = (\rho\theta)_k/\rho_k$ the horizontal part of the flux
divergence is discretized in the split form

$$
R^{\rho\theta,h}_k = -\tfrac12\Bigl(D_1(U^1\theta) + D_2(U^2\theta)
    + U^1 D_1\theta + U^2 D_2\theta + \theta\,(D_1 U^1 + D_2 U^2)\Bigr)_k \tag{28}
$$

(`DivRhoTrUpwind3New2Kernel!`). Together with the discrete continuity
equation this split form keeps a constant $\theta$ constant (the
split form reduces to $\theta$ times the mass flux divergence), and it is a
skew-symmetric (kinetic-energy-preserving type) discretization of the
advection of $\theta$. Tracers $\rho q_l$ use the plain flux form
$-(D_1(U^1 q_l) + D_2(U^2 q_l))$.

### 5.3 Vertical transport of rho theta and tracers

The vertical fluxes of $\rho\theta$ and $\rho q_l$ use the same mass flux
$\tfrac12 W_{k+1/2}$ as the continuity equation (consistency with mass) and
an upwind value of the specific quantity $c = \theta$ or $q_l$,

$$
\mathcal{F}_{k+1/2} = \tfrac14\Bigl((|W| + W)\,c^{L}_{k+1/2}
    + (W - |W|)\,c^{R}_{k+1/2}\Bigr), \qquad W = W_{k+1/2}. \tag{29}
$$

The interface values $c^{L}$, $c^{R}$ are reconstructed from the four layers
$k-1,k,k+1,k+2$ with a limited third-order upwind-biased reconstruction on
the non-uniform grid (`RecU4`). With the layer volumes
$J_k = \sum_K J^{(K)}_k$ and

$$
r = \frac{c_k - c_{k-1}}{c_{k+1} - c_k}, \qquad
  \kappa_R = \frac{J_k}{J_k+J_{k+1}}\,\frac{J_{k-1}+J_k}{J_{k-1}+J_k+J_{k+1}},\quad
  \kappa_L = \frac{J_k}{J_{k-1}+J_k}\,\frac{J_{k+1}}{J_{k-1}+J_k+J_{k+1}}, \tag{30}
$$

the left value is

$$
c^{L}_{k+1/2} = c_k + \max\bigl(0, \min(r,\, \kappa_L r + \kappa_R,\, 1)\bigr)\,(c_{k+1}-c_k), \tag{31}
$$

and $c^R$ is defined mirror-symmetrically. On a uniform grid
$\kappa_L = 1/6$, $\kappa_R = 1/3$, and without the limiter
$c^L = -\tfrac16 c_{k-1} + \tfrac56 c_k + \tfrac13 c_{k+1}$, the
third-order upwind-biased value. The limiter is of Koren type [Koren1993] (TVD, bounded by
$0$ and $1$ and by the first-order slope $r$). Near the boundaries the missing
layers are replaced by the boundary layer.

### 5.4 Momentum: rotational term

The vorticity is evaluated with the conservative curl (19):
the horizontal fluxes $\mathcal{J}^a\times\boldsymbol{v}$, $a=1,2$, are rotated to Cartesian
components before the derivative and the result is rotated back.
The horizontal derivatives are taken with $D$ separately at both vertical
nodes $K=1,2$ of the layer. The vertical derivative of a flux $f$ (for example
$\mathcal{J}^3_{\;3}v - \mathcal{J}^3_{\;2}w$) has two parts:

- inside the layer the two-point difference $\tfrac12(f^{(2)}_k - f^{(1)}_k)$,
- at each interface the jump $\tfrac12(f^{(1)}_{k+1} - f^{(2)}_k)$; each of
  the two adjacent half cells gets half of it (a central interface
  treatment, comparable to a DG lifting with the central flux).

The cross product with $\boldsymbol{v}$ is formed pointwise at the vertical nodes with the
interface values of $w$. The horizontal momentum tendency is

$$
\partial_t (u,v)_k \leftarrow
  -\frac{1}{M_{g,k}}\sum_{\text{copies}}\sum_{K=1}^2
  \Bigl[(J\boldsymbol\omega^{(K)})\times\boldsymbol{v}^{(K)}\Bigr]_{1,2} \tag{32}
$$

(`MomentumVectorInvariantCoriolisKernel!`). The vertical velocity gets
the corresponding contributions $(J\boldsymbol\omega^{(K)}\times\boldsymbol{v}^{(K)})_3$ of the
two half cells adjacent to an interface (Section 5.6).

The Coriolis force is added in the same kernel, in each half cell with its
own vertical velocity: `CoriolisFun` is evaluated at the node
$\boldsymbol X^{(K)}$ with $(u_k,v_k,w_{k\mp1/2})$ and its result is multiplied by
$J^{(K)}$, for $u$, $v$ and for the $w$ contribution of the half cell.
With
`CoriolisShallow` it returns $(f_sv,\,-f_su,\,0)$,
$f_s = 2\Omega\sin\varphi$, $\sin\varphi = z/|\boldsymbol X|$; with `FPlane`
$f_s = f_0 + \beta(y-y_0)$; with `CoriolisNonTraditional`

$$
\bigl(f_sv - f_cw,\; -f_su,\; f_cu\bigr),\qquad
  f_c = 2\Omega\cos\varphi . \tag{33}
$$

The non-traditional terms $-f_cw$ and $f_cu$ cancel in the energy within each
half cell, $u\,(-f_cw_{k\mp1/2}) + w_{k\mp1/2}\,f_cu = 0$, with the same
weights as the rotational term (the $w$ contributions are weighted with
$\rho_k$, Section 5.6); discretely the Coriolis contribution to the
kinetic energy vanishes to round-off. Together with the shallow metric this is
the NHNT equation set of [Hughes2026], with the deep metric their NHD
set.

### 5.5 Momentum: kinetic energy, pressure gradient and gravity

The kinetic energy at the two vertical nodes of a layer is

$$
K^{(1)}_k = \tfrac12\bigl(u_k^2 + v_k^2 + w_{k-1/2}^2\bigr),\qquad
  K^{(2)}_k = \tfrac12\bigl(u_k^2 + v_k^2 + w_{k+1/2}^2\bigr). \tag{34}
$$

Its gradient $J\nabla K = \sum_a\mathcal{J}^a\partial_{\xi^a}K$ is formed with $D$ in
the horizontal and with the same two-part vertical difference as above
(`GradFullKernel!`).

The horizontal pressure gradient uses the horizontal derivatives of the
layer pressure $p_k$:

$$
(J\nabla p)_b \approx \sum_{K=1}^2\sum_{a=1}^2 \mathcal{J}^{a,(K)}_{\;b}\,\partial_{\xi^a}p_k
  + \Bigl(\sum_{K}\mathcal{J}^{3,(K)}_{\;b}\Bigr)\,\partial_\zeta p,
  \qquad
  \partial_\zeta p \approx -g\,\rho_k\,
  \frac{J^{(1)}_k + J^{(2)}_k}{\mathcal{J}^{3,(1)}_{3,k} + \mathcal{J}^{3,(2)}_{3,k}},
  \qquad b = 1,2. \tag{35}
$$

Because $\mathcal{J}^3_{\;3} = J\,\partial\zeta/\partial z$, the quotient is
$\partial z/\partial\zeta$, so the vertical pressure derivative in the
terrain-following correction term is replaced by its hydrostatic value
$\partial_z p = -\rho g$. This is the treatment of the horizontal pressure
gradient over terrain used in MPAS [KlempSkamarock2021]. It needs no
vertical average of $\partial_\zeta p$ from the interfaces to the layer
centres, and it is exact for hydrostatic states. The horizontal tendency is
$-(J\nabla p)_{1,2}/\rho_k$, assembled with $M_{g,k}$. Without orography
$\mathcal{J}^3_{\;1} = \mathcal{J}^3_{\;2} = 0$ and the correction term vanishes; the first
term is then, by (14), the spherical gradient
$\bigl((a\cos\varphi)^{-1}\partial_\lambda p,\; a^{-1}\partial_\varphi p\bigr)$.

**Split form of the pressure gradient.**
With the option `--PressureGradient="Split"`
(kernel `GradFullSplitKernel!`, dry $\rho\theta$ only) the first term is replaced by the form that matches the
split transport of $\rho\theta$, as in the entropy stable formulation of
HOMME [TaylorSplit2025]. The derivative of the internal energy with
respect to $\rho\theta$ is $c_p\Pi$; multiplying the split transport by
$c_p\Pi$ and summing by parts gives the pressure term per unit mass

$$
\tfrac12 c_p\bigl(\theta\nabla\Pi + \nabla(\theta\Pi) - \Pi\nabla\theta\bigr)
  = \tfrac{c_p}2\Bigl(\theta\,\nabla\frac T\theta + \nabla T
    - \frac T\theta\,\nabla\theta\Bigr),\qquad
  T = \frac{p}{R_d\rho},\quad \Pi = \frac T\theta , \tag{36}
$$

with the horizontal SEM derivatives $D_a$ and $\sum_K\mathcal{J}^{a,(K)}_{\;b}$ as
in the transport. The Exner pressure is never evaluated: $p = R_d\rho\theta\Pi$
gives $\Pi = p/(R_d\rho\theta)$ without a fractional power. In the continuum
the expression equals $c_p\theta\nabla\Pi = \nabla p/\rho$; discretely the
horizontal pressure work and the horizontal $\rho\theta$ transport cancel to
round-off in the energy (relative $10^{-14}$), while the standard form
$\rho^{-1}D_ap$ leaves a truncation error ($2\cdot10^{-6}$ to $6\cdot10^{-10}$
for $n_{\text{Panel}} = 4$ to $16$). The terrain correction and the vertical
momentum equation are the same as in the standard form.

### 5.6 Vertical momentum

All contributions to the vertical velocity are formed in the two half cells
that touch the interface $k+\tfrac12$ (top half of layer $k$, bottom half of
layer $k+1$) and assembled with density-weighted masses:

$$
\partial_t w_{k+1/2} =
  \frac{\sum_{\text{copies}}\bigl(\rho_k R^{w,(2)}_k + \rho_{k+1} R^{w,(1)}_{k+1}
    + R^{p}_{k+1/2}\bigr)}
  {\sum_{\text{copies}}\bigl(\rho_k J^{(2)}_k + \rho_{k+1} J^{(1)}_{k+1}\bigr)}. \tag{37}
$$

Here $R^{w,(K)}$ contains the rotational term and the vertical gradient of
$K$ and, with the non-traditional Coriolis force, the term $f_cu$. The
vertical pressure gradient and gravity are

$$
R^p_{k+1/2} = -\tfrac12\bigl(p_{k+1} - p_k\bigr)
     \bigl(\mathcal{J}^{3,(2)}_{3,k} + \mathcal{J}^{3,(1)}_{3,k+1}\bigr)
     - g\bigl(\rho_k J^{(2)}_k + \rho_{k+1} J^{(1)}_{k+1}\bigr), \tag{38}
$$

i.e. $\bar\rho\,\partial_t w = -\partial_z p - \bar\rho g$ with the
volume-weighted interface density $\bar\rho$. Here $g$ is
`GravitationFun` evaluated at the interface; in the layer term of the
horizontal pressure gradient it is evaluated at the layer centre.
Because $\mathcal{J}^3_{\;3}$ is the radial component of $J\nabla\zeta$ in the
spherical frame, $p_{k+1}-p_k$ is differenced along the radial direction
also over orography. For a hydrostatic state
$(p_{k+1}-p_k)/\Delta z_{k+1/2} = -g\bar\rho_{k+1/2}$ is the discrete
balance that the right-hand side preserves exactly.

### 5.7 Discrete kinetic energy and vertical averaging

The definition of the discrete kinetic energy and the averaging to the
interfaces are chosen such that the vertical discretization conserves the
density-weighted kinetic energy of a column. We show this for the
two-dimensional system in the coordinates $(\xi,\zeta)$ written in covariant
($u_\xi,u_\zeta$) and contravariant ($u^\xi,u^\zeta$) components,

$$
\begin{aligned}
  J\partial_t\rho &= -\partial_\xi(J\rho u^\xi) - \partial_\zeta(J\rho u^\zeta), \qquad (39) \\
  J\rho\,\partial_t u_\xi &= -J\rho\,\partial_\xi K
     - J\rho u^\zeta\,\partial_\xi u_\zeta + J\rho u^\zeta\,\partial_\zeta u_\xi, \qquad (40) \\
  J\rho\,\partial_t u_\zeta &= -J\rho\,\partial_\zeta K
     + J\rho u^\xi\,\partial_\xi u_\zeta - J\rho u^\xi\,\partial_\zeta u_\xi, \qquad (41)
\end{aligned}
$$

with $K = \tfrac12(u_\xi u^\xi + u_\zeta u^\zeta)$; pressure gradient and
gravity are omitted. The covariant $u_k$ (and $u^\xi$) live in the layers,
$w_{k+1/2}$ (and $u^\zeta$) on the interfaces.

**Half-cell contravariant velocities.**
With the symmetric metric tensor $G_k$ of layer $k$, which maps covariant to
contravariant components, the velocities of the lower (L) and upper (R) half
cell are

$$
\begin{pmatrix} U^L_k\\ W^L_k\end{pmatrix}
    = G_k\begin{pmatrix} u_k\\ w_{k-1/2}\end{pmatrix},\qquad
  \begin{pmatrix} U^R_k\\ W^R_k\end{pmatrix}
    = G_k\begin{pmatrix} u_k\\ w_{k+1/2}\end{pmatrix}. \tag{42}
$$

The boundary values $w_{1/2}$ and $w_{N+1/2}$ are chosen such that
$W^L_1 = W^R_N = 0$. The layer and interface values are

$$
U_k = \tfrac12\bigl(U^L_k + U^R_k\bigr),\qquad
  W_{k+1/2} = \frac{J_k\rho_k W^R_k + J_{k+1}\rho_{k+1}W^L_{k+1}}
                   {J_k\rho_k + J_{k+1}\rho_{k+1}}, \tag{43}
$$

i.e. the contravariant vertical velocity at an interface is the
mass-weighted mean of the two adjacent half cells.

**Kinetic energy.**
The kinetic energy of a layer is the mean of the two half-cell energies,

$$
\begin{split}
  K_k &= \frac14\begin{pmatrix} u_k\\ w_{k-1/2}\end{pmatrix}^{\!T}
          G_k\begin{pmatrix} u_k\\ w_{k-1/2}\end{pmatrix}
      + \frac14\begin{pmatrix} u_k\\ w_{k+1/2}\end{pmatrix}^{\!T}
          G_k\begin{pmatrix} u_k\\ w_{k+1/2}\end{pmatrix}\\
      &= \frac14\bigl(u_k(U^L_k+U^R_k) + w_{k-1/2}W^L_k + w_{k+1/2}W^R_k\bigr).
\end{split} \tag{44}
$$

Because $G_k$ is symmetric,
$\dot K_k = \tfrac12\bigl(\dot u_k(U^L_k+U^R_k) + \dot w_{k-1/2}W^L_k + \dot w_{k+1/2}W^R_k\bigr)$.

**Semi-discrete equations.**
Let $D_x$ be the horizontal derivative, $D_{k+1/2}$ a discretization of
$\partial_\zeta u_\xi$ at the interface and
$\Omega_{k+1/2} = D_{k+1/2} - D_x w_{k+1/2}$ the discrete vorticity there.
Then

$$
\begin{aligned}
  J_k\dot\rho_k &= -D_x(J_k\rho_kU_k)
    - \tfrac12(J_k\rho_k + J_{k+1}\rho_{k+1})W_{k+1/2} \\
    &\quad + \tfrac12(J_{k-1}\rho_{k-1} + J_k\rho_k)W_{k-1/2}, \qquad (45) \\
  J_k\rho_k\dot u_k &= -J_k\rho_k D_xK_k
    + \tfrac12 J_k\rho_k\bigl(W_{k+1/2}\Omega_{k+1/2}
    + W_{k-1/2}\Omega_{k-1/2}\bigr), \qquad (46) \\
  \tfrac12(J_k\rho_k + J_{k+1}\rho_{k+1})\,\dot w_{k+1/2} &=
    -\tfrac12(J_k\rho_k + J_{k+1}\rho_{k+1})(K_{k+1} - K_k) \\
    &\quad - \tfrac12(J_k\rho_kU_k + J_{k+1}\rho_{k+1}U_{k+1})\,\Omega_{k+1/2}. \qquad (47)
\end{aligned}
$$

By (43), the vertical mass flux in (45) is
$\tfrac12(J_k\rho_kW^R_k + J_{k+1}\rho_{k+1}W^L_{k+1})$.

**Column energy.**
Summing $J_k\rho_kK_k$ over the column and using (43) and the
boundary conditions gives

$$
\mathcal K = \sum_{k=1}^N J_k\rho_kK_k
   = \sum_{k=1}^N J_k\rho_k u_kU_k
   + \sum_{k=1}^{N-1}\tfrac12(J_k\rho_k + J_{k+1}\rho_{k+1})\,w_{k+1/2}W_{k+1/2}, \tag{48}
$$

and therefore

$$
\begin{split}
  \dot{\mathcal K} = {}&\sum_{k=1}^N J_k\dot\rho_kK_k
   + \sum_{k=1}^N J_k\rho_k\dot u_kU_k\\
   &+ \sum_{k=1}^{N-1}\tfrac12(J_k\rho_k + J_{k+1}\rho_{k+1})\,\dot w_{k+1/2}W_{k+1/2}.
\end{split} \tag{49}
$$

Inserting (45)–(47), the vertical terms cancel pairwise:

- The vorticity terms cancel between the $u$ and the $w$ equation.
  This requires the same interface velocity $W_{k+1/2}$ in (46) as
  in the energy, and the mass-weighted layer velocity in (47).
- The vertical mass flux in the continuity equation times $K_k$ cancels,
  after summation by parts, against the vertical kinetic energy gradient
  in (47). This requires the same weight
  $\tfrac12(J_k\rho_k + J_{k+1}\rho_{k+1})$ in both places and the layer
  mean (44) for $K$.

What remains are boundary terms proportional to $W_{1/2}$ and $W_{N+1/2}$,
which vanish. The horizontal terms are handled by the SEM summation-by-parts
property of $D_x$ with the lumped GLL mass.

**Relation to the implementation.**
The code uses the orthonormal spherical components $(u,v,w)$. With respect
to an orthonormal basis covariant and contravariant components coincide, so
$G_k$ is the identity in these components, $K = \tfrac12(u^2+v^2+w^2)$, and
the metric enters only through the adjugate $\mathcal{J}$:

- The half-cell energies $K^{(1)}_k$, $K^{(2)}_k$ of
  `GradFullKernel!` are the two terms in (44), i.e.\
  $K_k = \tfrac12(K^{(1)}_k + K^{(2)}_k)$. $K_k$ itself is never formed.
  The horizontal gradient is applied to each half-cell energy with the
  metric of its node and weighted with the half-cell masses.
- The vertical mass flux of `Contra3` is exactly
  $\tfrac12(J_k\rho_kW^R_k + J_{k+1}\rho_{k+1}W^L_{k+1})$ with
  $J_kW^R_k = \mathcal{J}^{3,(2)}_k\cdot\boldsymbol{v}^{(2)}_k$, i.e. the averaging
  (43). The same flux is used for $\rho\theta$ and the tracers.
- The rotational term is formed pointwise in each half cell with the
  half-cell velocity $\boldsymbol{v}^{(K)}_k$. Since $(\boldsymbol\omega\times\boldsymbol{v})\cdot\boldsymbol{v} = 0$
  for any $\boldsymbol\omega$ (the shallow curvature correction is of the same form) and the
  $w$ contributions are weighted with $\rho_k$ and $\rho_{k+1}$ in the
  denominator of Section 5.6, it is energy neutral half cell by half
  cell. This takes the role of the averaged $W_{k+1/2}$ in (46).
- The vertical kinetic energy gradient in the $w$ equation is assembled
  half cell by half cell. In flat coordinates this gives

  $$
\tfrac12\rho_k\bigl(K^{(1)}_{k+1} - K^{(1)}_k\bigr)
    + \tfrac12\rho_{k+1}\bigl(K^{(2)}_{k+1} - K^{(2)}_k\bigr) \tag{50}
$$

  instead of $\tfrac12(\rho_k + \rho_{k+1})(K_{k+1} - K_k)$ in
  (47). The difference is
  $\tfrac14(\rho_k - \rho_{k+1}) \bigl[(K^{(1)}_{k+1} - K^{(1)}_k) - (K^{(2)}_{k+1} - K^{(2)}_k)\bigr]$.
  It vanishes for $\rho_k = \rho_{k+1}$ but not in general, so the column
  energy identity above holds for the code only up to this term.

### 5.8 Horizontal hyperviscosity

The dissipation $\boldsymbol D_{\boldsymbol{v}}$, $D_\theta$, $D_{q_l}$ is a biharmonic operator
in the horizontal, built in two passes with a DSS in between
(`HyperViscKernel!`, `HyperViscKoeffKernel!`).

1. First pass: per layer the vertical vorticity
  $\zeta_v = \boldsymbol k\cdot\nabla\times\boldsymbol{v}$ and the divergence
  $\delta = \nabla_h\cdot\boldsymbol{v}$ are computed in strong form and divided by
  $J^{(1)}+J^{(2)}$. Then the weak forms of the rotational and the divergent
  part of the vector Laplacian,

  $$
\boldsymbol L_{\text{rot}} = -\boldsymbol k\times\nabla^W \zeta_v,\qquad
    \boldsymbol L_{\text{div}} = \nabla^W\delta,\qquad
    L_\theta = \nabla^W\cdot\nabla\theta, \tag{51}
$$

  are formed with $D^W$ (so that $\boldsymbol L_{\text{rot}}+\boldsymbol L_{\text{div}} = \nabla(\nabla\cdot\boldsymbol{v}) - \nabla\times(\nabla\times\boldsymbol{v})$ is the vector
  Laplacian). They are stored separately and assembled (DSS + exchange).
2. Second pass: the same operators are applied to the stored fields:

  $$
\boldsymbol D_{\boldsymbol{v}} = -\nu_{\text{curl}}\,\Delta_{\text{rot}}\boldsymbol L_{\text{rot}}
       - \nu_{\text{grad}}\,\Delta_{\text{div}}\boldsymbol L_{\text{div}},\qquad
    D_\theta = -\nu_{\theta}\,\nabla^W\cdot\bigl(\rho\,\nabla L_\theta\bigr). \tag{52}
$$

The coefficients are `HyperDCurl`, `HyperDGrad` and
`HyperDDiv`. Optionally there is a hyperviscosity for $w$
(`HyperDDivW`), a biharmonic of $w$ on the interfaces with the
interface metric $\mathcal{J}^{(2)}_k + \mathcal{J}^{(1)}_{k+1}$. Tracers (and the TKE) are
treated like $\theta$ with the specific value $q_l = (\rho q_l)/\rho$ in the
first pass and the density weight in the second pass; with
`HorLimit=true` the tracer hyperviscosity is switched off.

All hyperviscosity operators use only the horizontal entries
$\mathcal{J}^a_{\;b}$, $a,b = 1,2$ of the adjugate, summed over the two vertical
nodes; they act along the coordinate surfaces $\zeta = \text{const}$ and
contain no terrain correction. In the spherical frame the curl and the
divergence of the first pass are, by (14), the spherical
vertical vorticity and horizontal divergence on these surfaces. The default scaling
with the resolution (`HyperViscosity`) is

$$
\nu = 5\cdot 10^{14}\,\mathrm{m^4\,s^{-1}}\;
  \Bigl(\frac{32}{n_{\text{Panel}}}\Bigr)^{\log_2 10}
  \approx 5\cdot 10^{14}\,\Bigl(\frac{32}{n_{\text{Panel}}}\Bigr)^{3.32}, \tag{53}
$$

comparable to the $\Delta x^{3.2}$ scaling that is used in HOMME.

### 5.9 Optional horizontal limiter for tracers

With `HorLimit=true` the horizontal tracer transport uses an
optimization-based limiter in the spirit of Guba et al. (2014)
(`LimitKernel!`, `DivRhoTrUpwind3LimKernel!`):

1. the minimum and maximum of $q$ in every element are computed, and the
  bounds $q_{\min}, q_{\max}$ are taken over the element and its neighbours
  (a stencil of 13 elements);
2. a provisional explicit step $\tilde\rho$, $\widetilde{\rho q}$ with the
  current stage time step $\Delta\tau$ (`dtauStage`) is formed from the
  element-local horizontal flux divergences (before DSS, without the
  vertical flux);
3. a constant shift $\lambda$ is determined per element and layer so that
  $q^\star = \operatorname{median}(q_{\min}, \widetilde{\rho q}/\tilde\rho + \lambda, q_{\max})$
  conserves the tracer mass in the element,
  $\sum \omega_i\omega_j J\,(q^\star\tilde\rho - \widetilde{\rho q}) = 0$.
  The scalar equation is solved by a secant iteration (at most 8 steps).

The limited horizontal tendency of the element,
$J\,(q^\star\tilde\rho - \rho q)/\Delta\tau$, is then assembled with DSS like
every other term. The vertical part uses the reconstruction of
Section 5.3 and is not limited. In this mode $\rho\theta$ is
transported in plain flux form (`DivRhoThUpwind3Kernel!`) instead of
the split form.

## 6 Cartesian velocity frame (`--VelocityFrame=Cartesian`)

As an alternative to the spherical components $(u,v,w)$ the horizontal
velocity can be stored with its Cartesian components. The prognostic
velocity is

$$
\boldsymbol{v} = \boldsymbol u_h + w\,\hat{\boldsymbol n},\qquad
  \boldsymbol u_h = (u_x,u_y,u_z)^T,\quad \hat{\boldsymbol n}\cdot\boldsymbol u_h = 0,\qquad
  \hat{\boldsymbol n} = \frac{\boldsymbol X}{|\boldsymbol X|}, \tag{54}
$$

with $\boldsymbol u_h$ at the positions `uPos`, `vPos` and
`uzPos` $=6$ (`NumV` $\ge 6$, the tracers follow `NumV`) and
$w$ at the interfaces as before. The rotation matrix $R$ of
Section 3.3 is never formed: the only quantity that depends
on the sphere besides $\boldsymbol X$ is $\hat{\boldsymbol n}$, which costs one square root.
No trigonometric function is evaluated at run time; $\sin$ and $\cos$ only
appear in the conversion of the initial state and of the output
(`ToVelocityFrame!`, `FromVelocityFrame`). The spherical frame
remains the default (`Models.SphericalFrame`); the new option is
`Models.CartesianFrame(uzPos)`.

**Metric.**
`dXdxI` gets a fourth column. Columns $1$–$3$ hold the unrotated
Cartesian rows $\mathcal{J}^a_C$ and column $4$ the normal part

$$
\texttt{dXdxI}[a,1{:}3] = \mathcal{J}^a_C,\qquad
  \texttt{dXdxI}[a,4] = \mathcal{J}^a_C\cdot\hat{\boldsymbol n}, \tag{55}
$$

(`NormalColumnKernel!` instead of `ApplyRotateKernel!`), so that the
contravariant velocity is
$\mathcal{J}^a_C\cdot\boldsymbol{v} = \mathcal{J}^a_C\cdot\boldsymbol u_h + (\mathcal{J}^a_C\cdot\hat{\boldsymbol n})\,w$.
The shallow metric is the Cartesian form of $L\to SL$: with $s = a/r$,

$$
\mathcal{J}^a_C \to s\,\mathcal{J}^a_C + (s^2-s)\,(\hat{\boldsymbol n}\cdot\mathcal{J}^a_C)\,\hat{\boldsymbol n},\qquad
  \hat{\boldsymbol n}\cdot\mathcal{J}^a_C \to s^2\,\hat{\boldsymbol n}\cdot\mathcal{J}^a_C,\qquad
  J\to s^2 J \tag{56}
$$

(`ShallowMetricCartKernel!`), i.e. the tangential part is scaled by
$s$ and the normal part by $s^2$, exactly as $s = (a/r,a/r,(a/r)^2)$ in the
spherical frame. The node normals are stored once per DSS node in
`Metric.nR`. The lower boundary condition $\mathcal{J}^3\cdot\boldsymbol{v} = 0$ gives

$$
w_{1/2} = -\frac{\mathcal{J}^3_C\cdot\boldsymbol u_h}{\mathcal{J}^3_C\cdot\hat{\boldsymbol n}} \tag{57}
$$

(`BottomW`).

**Tangential and normal part of the rows.**
After $J$, the shallow scaling and the surface normals, which use the full
rows, are computed, columns $1$–$3$ keep only the tangential part
(`TangentialRows!`),

$$
\texttt{dXdxI}[a,1{:}3] = \mathcal{J}^a_C - (\mathcal{J}^a_C\cdot\hat{\boldsymbol n})\,\hat{\boldsymbol n},
  \qquad \texttt{dXdxI}[a,4] = \mathcal{J}^a_C\cdot\hat{\boldsymbol n}. \tag{58}
$$

Since $\boldsymbol u_h$ is tangential, the contravariant velocity and the lower
boundary condition are unchanged. A gradient
$\sum_a\mathcal{J}^a_C\partial_{\xi^a}\phi$ now contributes its tangential part to
$\boldsymbol u_h$ directly. With the full rows, the large normal part of
$\mathcal{J}^3_C$ (gravity and the vertical pressure and kinetic energy gradients,
$\approx 10\,$m s$^{-2}$) entered the tendency of $\boldsymbol u_h$ and was removed
again by the projection below. In Float32 this cancellation left a residual
of relative size $10^{-7}$ with the same sign in every step, $\hat{\boldsymbol n}\cdot\boldsymbol u_h$
grew by $10^{-2}\,$m s$^{-1}$ per hour and the baroclinic wave became unstable
after one day. With the split rows the growth is $4.5\cdot10^{-5}\,$m s$^{-1}$
per hour and the Float32 solution after two days deviates from the Float64 one as
much as in the spherical frame. Only the cross products of the curl need the
full rows, which are rebuilt as $\texttt{dXdxI}[a,1{:}3] + \texttt{dXdxI}[a,4]\,\hat{\boldsymbol n}$.

**Operators.**
Divergence and transport (`DivRhoTrUpwind3New2CartKernel!`), the
pressure gradient in standard and split form, kinetic energy and gravity
(`GradFullCartKernel!`) and the hyperviscosity
(`HyperViscCartKernel!`, `HyperViscKoeffCartKernel!`) are the
operators (10) with $\mathcal{J}^a_C$. A gradient contributes the full
Cartesian vector $\sum_a\mathcal{J}^a_C\partial_{\xi^a}\phi$ to $\boldsymbol u_h$ and its
normal component (column $4$) to $w$. In the hyperviscosity the vorticity
is $\zeta = \hat{\boldsymbol n}\cdot\sum_a\partial_{\xi^a}(\mathcal{J}^a_C\times\boldsymbol u_h)$,
with the normal product formed before the derivative, and the cache of the
first pass holds the three Cartesian components of $\hat{\boldsymbol n}\times\nabla\zeta$
and of $\nabla\delta$ (seven instead of five slots).

**Rotational term.**
At each vertical node $K$ of a layer the full velocity
$\boldsymbol{v}^{(K)} = \boldsymbol u_h + w^{(K)}\hat{\boldsymbol n}$ is formed with the interface value
$w^{(K)} = w_{k\mp1/2}$. The curl is the Cartesian curl (10),

$$
J\boldsymbol\omega^{(K)} = \sum_{a=1}^2 D_a\bigl(\mathcal{J}^a_C\times\boldsymbol{v}\bigr)^{(K)}
   + D_\zeta\bigl(\mathcal{J}^3_C\times\boldsymbol{v}\bigr), \tag{59}
$$

with the same two-point difference and interface jumps as in
Section 5.4. No `CurlFrame` rotation is needed, since the
Cartesian basis is constant. The tendency
$\boldsymbol A^{(K)} = -J\boldsymbol\omega^{(K)}\times\boldsymbol{v}^{(K)}$ is added to $\boldsymbol u_h$ and
$\hat{\boldsymbol n}\cdot\boldsymbol A^{(K)}$ to the $w$ contribution of the half cell. With
the shallow metric the curvature terms are removed as in the spherical
frame,

$$
\partial_t\boldsymbol u_h \leftarrow -\frac{1}{a}\bigl(J^{(1)}w^{(1)} + J^{(2)}w^{(2)}\bigr)\boldsymbol u_h,
  \qquad
  \partial_t w^{(K)} \leftarrow \frac{1}{a}J^{(K)}|\boldsymbol u_h|^2 . \tag{60}
$$

The Coriolis functions have Cartesian versions that return the
tangential vector and the normal component:
`CoriolisShallow` $-2\Omega\hat n_z\,\hat{\boldsymbol n}\times\boldsymbol u_h$, `FPlane`
the same with $f_0+\beta(y-y_0)$, and `CoriolisNonTraditional`
$\boldsymbol A = -2\Omega\,\hat{\boldsymbol z}\times(\boldsymbol u_h + w\hat{\boldsymbol n})$ with
$\hat{\boldsymbol n}\cdot\boldsymbol A$ for $w$.

**Tangential projection.**
The kernels add unprojected Cartesian vectors to the tendency of $\boldsymbol u_h$.
At the end of `Fcn!`, after DSS and the division by the mass,
`ProjectTangential!` applies once per node

$$
F_{\boldsymbol u_h} \leftarrow F_{\boldsymbol u_h} - (\hat{\boldsymbol n}\cdot F_{\boldsymbol u_h})\,\hat{\boldsymbol n}, \tag{61}
$$

so that $\hat{\boldsymbol n}\cdot\boldsymbol u_h = 0$ is preserved by every stage. Since the
normal parts already enter the $w$ equation through column $4$ and
$\hat{\boldsymbol n}\cdot\boldsymbol A$, the projection only removes the redundant normal
component; it does not change the dynamics. After every time step the state
itself is projected in the same way (`Integration.ProjectState!`). The
update $\boldsymbol u_h + \sum_i m_i\boldsymbol k_i$ rounds the three Cartesian components
separately, and in a nearly steady flow the rounding has the same sign step by
step: in Float32 $\hat{\boldsymbol n}\cdot\boldsymbol u_h$ grew linearly to $2\cdot10^{-3}\,$m s$^{-1}$
in two days, acted like an additional $w$ in the curl term and gave a
deviation of $0.16$ in $w$ from the Float64 solution. With the projection
$\hat{\boldsymbol n}\cdot\boldsymbol u_h$ stays at the rounding level ($2\cdot10^{-6}\,$m s$^{-1}$) and
the Float32 solution of the Cartesian frame is as accurate as the one of the
spherical frame.

**Time integration.**
The vertical Jacobian (Section 7.2) depends only on $\Delta z$ and
is unchanged; $u_z$ is treated explicitly like $u$ and $v$.

**Scope and verification.**
The option is restricted to the dry core
(`CGSEM.CheckVelocityFrame`): surface fluxes, momentum diffusion,
forcing, turbulence, EDMF, the horizontal limiter and the energy states still
read spherical $u$, $v$ and are rejected. Compared with the spherical frame
on the baroclinic wave ($n_{\text{Panel}} = 6$, $p = 3$, $12$ levels), every
term agrees to $10^{-13}$ relative (deep and shallow metric, with
orography, standard and split pressure gradient, non-traditional Coriolis,
hyperviscosity). After one day of integration $\rho$ and $\rho\theta$ agree
to $2\cdot10^{-13}$ and $u$, $v$ to $10^{-11}$; after ten days the
difference has grown to $5\cdot10^{-7}$ by the baroclinic instability.

## 7 Time integration

### 7.1 Rosenbrock W-methods

The semi-discrete system $U' = F(U)$ is integrated with an $s$-stage
Rosenbrock-type method in the transformed form [HairerWanner]
(`Integration/Rosenbrock.jl`)

$$
\begin{aligned}
  Y_i &= U^n + \sum_{j<i} a_{ij}\,k_j, \qquad (62) \\
  \Bigl(\frac{1}{\gamma\Delta t}I - \mathbf J\Bigr)k_i &= F(Y_i)
    + \frac{1}{\Delta t}\sum_{j<i} c_{ij}\,k_j, \qquad i = 1,\dots,s, \qquad (63) \\
  U^{n+1} &= U^n + \sum_{i} m_i\,k_i, \qquad (64)
\end{aligned}
$$

with $A = \alpha\Gamma^{-1}$, $C = -\Gamma^{-1}$ (off-diagonal part) and
$m = \Gamma^{-T}b$. The Jacobian approximation $\mathbf J$ is computed once
per step from the state $U^n$ in the first stage (W-method,
`CallJac!`) and contains only the vertical
coupling terms (HEVI), so the method is explicit in the horizontal and
linearly implicit in the vertical. The table `SSP-Knoth` has
$s=3$, $\gamma=1$,

$$
\alpha = \begin{pmatrix} 0&0&0\\ 1&0&0\\ \tfrac14&\tfrac14&0\end{pmatrix},\quad
  b = \begin{pmatrix}\tfrac16\\ \tfrac16\\ \tfrac23\end{pmatrix},\quad
  \Gamma = \begin{pmatrix} 1&0&0\\ 0&1&0\\ -\tfrac34&-\tfrac34&1\end{pmatrix}; \tag{65}
$$

its explicit part ($\alpha$, $b$) is the three-stage SSP Runge–Kutta method
of Shu and Osher [ShuOsher1988]. Its stability function tends to
$-1/6$ for $z\to-\infty$. The table `SSP-KnothL` [Artiano2026]
has the same $\alpha$ and $b$, $\gamma = 1-1/\sqrt2$ and

$$
\Gamma = \begin{pmatrix} \gamma&0&0\\ \tfrac3{\sqrt2}-3&\gamma&0\\
    \tfrac{3\sqrt2-4}8&-\tfrac14&\gamma\end{pmatrix}; \tag{66}
$$

it is a W-method of order 2 and L-stable. Further tables (`SSPKnoth12L`, `ROS2W`,
`ROSRK3`, `ARS343`, …) and explicit Runge–Kutta methods are
available. Optionally (`CompensatedSum`) the final update
$U^n + \sum_i m_i k_i$ is accumulated with Kahan compensated summation.

### 7.2 Vertical Jacobian and Schur complement

Per column and horizontal node the Jacobian couples $w$ with $\rho$ and
$\rho\theta$ (`JacAcousticKernel!`). It is built from the geometric
layer thicknesses $\Delta z_k$ (radial distance of the two vertical nodes,
`Metric.dz`) and the constant $g$, not from the metric terms of the
spatial operator; the terrain slope and the variation of $g$ with height are
therefore not contained in $\mathbf J$. For the total-energy variant the
role of $\rho\theta$ in the interface value is taken by the total enthalpy
$\rho e + p$, and $\partial_\rho p \neq 0$.
With the interface values
$\rho_F = (\rho_k\Delta z_k + \rho_{k+1}\Delta z_{k+1})/(\Delta z_k+\Delta z_{k+1})$
(and analogously $(\rho\theta)_F$) and $\Delta z_{k+1/2} = \tfrac12(\Delta z_k + \Delta z_{k+1})$, the non-zero blocks are

$$
\begin{aligned}
  \frac{\partial F_{\rho,k}}{\partial w_{k\pm1/2}} &= \mp\frac{\rho_{F,k\pm1/2}}{\Delta z_k},
  &
  \frac{\partial F_{\rho\theta,k}}{\partial w_{k\pm1/2}} &=
    \mp\frac{(\rho\theta)_{F,k\pm1/2}}{\Delta z_k}, \qquad (67) \\
  \frac{\partial F_{w,k+1/2}}{\partial (\rho\theta)_{k+1}} &=
    -\frac{\partial_{\rho\theta}p_{k+1}}{\rho_F\,\Delta z_{k+1/2}},
  &
  \frac{\partial F_{w,k+1/2}}{\partial (\rho\theta)_{k}} &=
    \frac{\partial_{\rho\theta}p_{k}}{\rho_F\,\Delta z_{k+1/2}}, \qquad (68) \\
  \frac{\partial F_{w,k+1/2}}{\partial \rho_{k+1}} &=
    -\frac{g\,\Delta z_{k+1}}{\rho_F(\Delta z_k+\Delta z_{k+1})},
  &
  \frac{\partial F_{w,k+1/2}}{\partial \rho_{k}} &=
    -\frac{g\,\Delta z_{k}}{\rho_F(\Delta z_k+\Delta z_{k+1})} \qquad (69)
\end{aligned}
$$

(plus $\partial_\rho p$ terms, which vanish for the dry $\rho\theta$
formulation). The linear system

$$
\begin{pmatrix}
    \tfrac1{\gamma\Delta t}I & 0 & -J_{\rho w}\\
    0 & \tfrac1{\gamma\Delta t}I & -J_{\rho\theta w}\\
    -J_{w\rho} & -J_{w\rho\theta} & \tfrac1{\gamma\Delta t}I
  \end{pmatrix}
  \begin{pmatrix} k_\rho\\ k_{\rho\theta}\\ k_w\end{pmatrix}
  =
  \begin{pmatrix} r_\rho\\ r_{\rho\theta}\\ r_w\end{pmatrix} \tag{70}
$$

is reduced to the Schur complement for $w$,

$$
\Bigl(\frac{1}{(\gamma\Delta t)^2}I - J_{w\rho}J_{\rho w}
    - J_{w\rho\theta}J_{\rho\theta w}\Bigr) k_w
  = \frac{1}{\gamma\Delta t} r_w + J_{w\rho} r_\rho + J_{w\rho\theta} r_{\rho\theta}, \tag{71}
$$

which is tridiagonal ($J_{\rho w}$, $J_{\rho\theta w}$ are lower and
$J_{w\rho}$, $J_{w\rho\theta}$ upper bidiagonal). It is solved per column
with the Thomas algorithm (`SchurSolveTriKernel!`), followed by the back
substitution $k_\rho = \gamma\Delta t\,(r_\rho + J_{\rho w}k_w)$,
$k_{\rho\theta} = \gamma\Delta t\,(r_{\rho\theta} + J_{\rho\theta w}k_w)$
and $k_{u,v,q} = \gamma\Delta t\,r_{u,v,q}$. The tridiagonal matrix is
assembled in the first solve of a step (flag `CompTri`) and reused in the
later stages, since $\gamma$ is the same for all stages. Optionally the vertical diffusion of the scalars
is also treated implicitly by an extra tridiagonal solve per tracer
(`SchurSolveTriDiffKernel!`).

## 8 Comparison with HOMME and ClimaCore

All three dynamical cores use the same horizontal discretization: continuous
Galerkin spectral elements on the cubed sphere with GLL collocation, lumped
mass, DSS and a vector-invariant momentum equation; they differ in the vertical
coordinate, the vertical discretization, the thermodynamic variable and the
time integration. Table 2 summarizes the differences. The
entries for HOMME and ClimaCore are based on the cited literature and
describe the standard configurations.

**Table 2:** Comparison of the three spectral element dynamical cores.

|  | CGDycore.jl (CGSEM) | HOMME / HOMME-NH (E3SM) | ClimaCore / ClimaAtmos (CliMA) |
|---|---|---|---|
| Horizontal | CG-SEM, GLL, $p=3$–$4$, DSS | CG-SEM, GLL, $n_p=4$ ($p=3$), DSS | CG-SEM, GLL, typically $p=3$, DSS |
| Vertical coordinate | terrain-following height | floating Lagrangian mass coordinate, remap to hybrid $\eta$ | terrain-following height |
| Vertical discretization | lowest-order staggered, $w$ on interfaces, metric from two-point layers | mimetic finite differences: $\Delta\pi$, $\Theta$ at mid levels, $w$, $\phi$ at interfaces | staggered finite differences (centers/faces) |
| Prognostic variables | $\rho$, $u,v,w$ (orthonormal spherical components), $\rho\theta$, $\rho q$ | $\boldsymbol u_h$, $w$, $\phi$, $\Delta\pi$, $\Theta=\theta_v\Delta\pi$, $q\,\Delta\pi$ | $\rho$, covariant $\boldsymbol u_h$, covariant $u_3$, $\rho e_{tot}$, $\rho q$ |
| Thermodynamics | potential temperature (energy variants) | virtual potential temperature | total energy |
| Vertical transport | limited third-order upwind (Koren type) | implicit through the Lagrangian coordinate and conservative remap | first/third-order upwind, optional flux-corrected transport |
| Dissipation | biharmonic, separate $\nu$ for curl/div/scalars | biharmonic (Laplacian or tensor), subcycled, sponge layer | biharmonic hyperdiffusion |
| Tracer limiter | optional optimization-based limiter | optimization-based limiter, or semi-Lagrangian transport | limiters / flux correction |
| Time integration | Rosenbrock W-methods (HEVI) | IMEX Runge–Kutta (HEVI) | IMEX additive Runge–Kutta (e.g. ARS343) |
| Implicit part | linear vertical acoustics and buoyancy, Schur complement for $w$, frozen Jacobian | vertical acoustics in $(w,\phi)$, Newton iteration | vertical acoustics (and optionally diffusion), Newton with approximate Jacobian |

**HOMME.**
HOMME (the High-Order Methods Modeling Environment, dynamical core of CAM-SE
and E3SM) introduced the CG-SEM on the cubed sphere for atmospheric models
[TaylorFournier2010, Dennis2012]. The nonhydrostatic version
[Taylor2020] uses a floating Lagrangian vertical coordinate in which the
vertical mass flux vanishes; the layers are remapped conservatively to a
reference hybrid coordinate after a number of steps. The vertical
discretization is mimetic and the semi-discrete scheme conserves the total
energy. Time integration is by IMEX Runge–Kutta methods [Steyer2019]
in which only the vertical acoustic terms are implicit; the resulting
nonlinear system for $w$ and $\phi$ is solved by Newton's method with
tridiagonal linear systems. Hyperviscosity and the limiter
[Guba2014a, Guba2014b] are applied in separate steps.

**ClimaCore.**
ClimaCore.jl provides the operators of the CliMA atmosphere model
[Yatunin2025]. It combines the horizontal SEM with a staggered
finite-difference discretization in a terrain-following height coordinate,
like CGDycore. The momentum is stored in covariant components (horizontal
on the cell centers, vertical on the faces) and the thermodynamic variable
is the total energy. The IMEX additive Runge–Kutta methods treat the vertical
acoustic terms implicitly with a Newton iteration and an approximate Jacobian.

**Specific features of CGDycore.**
Compared to the two other cores, the CGSEM scheme

- stores the velocity in orthonormal spherical components (east, north,
  radial) instead of covariant components, forms the Jacobian of the element
  mapping in the local spherical frame and evaluates all operators in the
  conservative (adjugate) form (15) at the two vertical
  nodes of each layer, so that the curvature terms of the sphere are
  contained in the metric and need no separate treatment;
- assembles the vertical velocity with density-weighted half-cell masses,
  which gives a $\bar\rho\,\partial_t w$ form of the vertical momentum
  equation and an exact discrete hydrostatic balance;
- uses a split (skew-symmetric) form for the horizontal advection of
  $\theta$ and a limited third-order upwind scheme in the vertical;
- uses linearly implicit Rosenbrock W-methods instead of IMEX Runge–Kutta
  methods: no Newton iteration is needed, the Jacobian is frozen over the
  step and the linear algebra reduces to one tridiagonal solve per column
  and stage;
- computes the vertical pressure derivative in the terrain-following
  correction of the horizontal pressure gradient from the hydrostatic
  relation.

## 9 Remarks

- The approximation $\partial_z p \approx -\rho g$ in the terrain term of
  the horizontal pressure gradient (as in MPAS [KlempSkamarock2021]) is
  exact for hydrostatic states. It neglects the vertical derivative of the
  nonhydrostatic pressure perturbation.
- The kinetic energy and the interface velocity are defined as in
  Section 5.7. The vertical kinetic energy gradient in the $w$
  equation is assembled half cell by half cell rather than from the layer
  mean $K_k$. This differs from the energy-conserving form by a term
  proportional to $\rho_k - \rho_{k+1}$.
- The Jacobian contains neither vertical advection nor the
  $u,v$-dependence of the lower boundary condition. Its accuracy influences
  only the stability and the order of the W-method, not the consistency of
  the scheme.
- Moving the rotation $R$ into the metric is exact for all fluxes,
  divergences and gradients, so the spherical scheme is node by node the
  rotated Cartesian scheme for these terms. For the curl the derivative of
  $R$ is not contained in the metric (16). The horizontal
  derivatives in the curl therefore act on the Cartesian components
  (19); with the shallow metric the remaining curvature
  terms $(uw,\,vw,\,-(u^2+v^2))/a$ are subtracted. Both are energy neutral
  and regular at the poles.
- The interface treatment of the vertical derivatives in the rotational
  term and in the kinetic energy gradient distributes the jump equally to the
  two adjacent half cells (central treatment). The vertical upwinding only
  acts on the scalar transport.

## References

- **[TaylorFournier2010]** M. A. Taylor and A. Fournier, A compatible and conservative spectral element method on unstructured grids, *J. Comput. Phys.* 229 (2010) 5879–5895.
- **[Dennis2012]** J. M. Dennis, J. Edwards, K. J. Evans, O. Guba, P. H. Lauritzen, A. A. Mirin, A. St-Cyr, M. A. Taylor, P. H. Worley, CAM-SE: A scalable spectral element dynamical core for the Community Atmosphere Model, *Int. J. High Perform. Comput. Appl.* 26 (2012) 74–89.
- **[Taylor2020]** M. A. Taylor, O. Guba, A. Steyer, P. A. Ullrich, D. M. Hall, C. Eldred, An energy consistent discretization of the nonhydrostatic equations in primitive variables, *J. Adv. Model. Earth Syst.* 12 (2020) e2019MS001783.
- **[Steyer2019]** A. Steyer, C. J. Vogl, M. Taylor, O. Guba, Efficient IMEX Runge–Kutta methods for nonhydrostatic dynamics, arXiv:1906.07219 (2019).
- **[Guba2014a]** O. Guba, M. A. Taylor, P. A. Ullrich, J. R. Overfelt, M. N. Levy, The spectral element method (SEM) on variable-resolution grids: evaluating grid sensitivity and resolution-aware numerical viscosity, *Geosci. Model Dev.* 7 (2014) 2803–2816.
- **[Guba2014b]** O. Guba, M. Taylor, A. St-Cyr, Optimization-based limiters for the spectral element method, *J. Comput. Phys.* 267 (2014) 176–195.
- **[Yatunin2025]** D. Yatunin, S. Byrne, C. Kawczynski, S. Kandala, G. Bozzola, A. Sridhar, Z. Shen, A. Jaruga, J. Solomon, C. Kao, et al., The CliMA atmosphere dynamical core: concepts, numerics, and scaling, preprint (2025).
- **[KlempSkamarock2021]** J. B. Klemp and W. C. Skamarock, Adapting the MPAS dynamical core for applications extending into the thermosphere, *J. Adv. Model. Earth Syst.* 13 (2021) e2021MS002499.
- **[Koren1993]** B. Koren, A robust upwind discretization method for advection, diffusion and source terms, in: *Numerical Methods for Advection–Diffusion Problems*, Vieweg, Braunschweig, 1993, pp. 117–138.
- **[HairerWanner]** E. Hairer and G. Wanner, *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*, 2nd ed., Springer, 1996.
- **[Hughes2026]** O. K. Hughes, O. Guba, M. A. Taylor, C. Jablonowski, The deep atmosphere extension of the non-hydrostatic HOMME dynamical core, *EGUsphere* preprint, doi:10.5194/egusphere-2026-4425 (2026).
- **[TaylorSplit2025]** M. A. Taylor, M. Waruszewski, C. Eldred, O. Guba, An entropy stable formulation of E3SM's spectral element atmospheric dynamical core, poster, E3SM All-Hands Meeting, August 2025.
- **[Artiano2026]** M. Artiano, A. M. Rueda-Ram\'irez, O. Knoth, H. Ranocha, G. J. Gassner, GPU-accelerated nested Schur-complement solvers for atmospheric flow simulations using horizontally-explicit, vertically-implicit schemes and discontinuous Galerkin methods, arXiv:2610.15775 (2026).
- **[ShuOsher1988]** C.-W. Shu and S. Osher, Efficient implementation of essentially non-oscillatory shock-capturing schemes, *J. Comput. Phys.* 77 (1988) 439–471.

