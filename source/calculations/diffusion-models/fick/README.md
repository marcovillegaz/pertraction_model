# fick/

Assembles the effective Fick diffusivity matrix from the Maxwell-Stefan
$[B]$ and $[\Gamma]$ matrices. Math lives in
[`docs/theory/diffusion-models/fick.md`](../../../../docs/theory/diffusion-models/fick.md).

- `computeFickDiffusivity.m` — entry point; calls `mutualDiffusion`, `Bmatrix`, `thermodynamicsFactors` and solves $[\text{Đ}] = [B]^{-1}[\Gamma]$
