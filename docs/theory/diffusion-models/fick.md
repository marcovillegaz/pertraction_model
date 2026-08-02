# Effective Fick diffusivity

Governs: `source/calculations/diffusion-models/fick/computeFickDiffusivity.m`.

Final assembly step of the diffusion chain (see the pipeline diagram in
[`../README.md`](../README.md)). Converts the generalized Maxwell-Stefan
formulation into the effective Fick diffusivity matrix used in the mass
transfer equation:

$$
[\text{Đ}] = [B]^{-1}[\Gamma]
$$

where:

- $[B]$ is the inverse drag matrix from mutual diffusion coefficients — see [`maxwell-stefan.md`](maxwell-stefan.md) § 2
- $[\Gamma]$ is the thermodynamic factor matrix — see [`maxwell-stefan.md`](maxwell-stefan.md) § 3
- $[\text{Đ}]$ is $(m-1)\times(m-1)$, units cm²/s

Upstream, $[B]$ requires the full mutual diffusion matrix, assembled from:

- non-polymer pairs — [`dilute-correlations.md`](dilute-correlations.md)
- polymer pairs — [`maxwell-stefan.md`](maxwell-stefan.md) § 1, which itself needs self-diffusion coefficients from [`free-volume-theory.md`](free-volume-theory.md)

**Reference:** Taylor & Krishna (1993). See [`../references.md`](../references.md).
