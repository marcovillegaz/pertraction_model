# maxwellStefan/

Polymer mutual diffusion, the inverse-drag matrix, and the thermodynamic
factor matrix. Math lives in
[`docs/theory/diffusion-models/maxwell-stefan.md`](../../../../docs/theory/diffusion-models/maxwell-stefan.md).

- `mutualDiffusion.m` — assembles the full Maxwell-Stefan diffusion matrix (non-polymer pairs + `kubaczka.m` polymer pairs)
- `thermodynamicsFactors.m` — computes $[\Gamma]$ via finite-difference Jacobian of $\ln\gamma$ (UNIFAC — see [`docs/theory/thermodynamic-models/unifac.md`](../../../../docs/theory/thermodynamic-models/unifac.md))
- `Bmatrix.m` — computes the inverse-drag matrix $[B]$
