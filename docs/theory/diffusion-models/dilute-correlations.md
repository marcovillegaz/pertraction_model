# Dilute and mutual diffusion correlations (non-polymer pairs)

Governs: `source/calculations/diffusion-models/siddiqiLucas.m`,
`source/calculations/diffusion-models/kooijmanTaylor.m`.

These two steps produce the non-polymer block of the Maxwell-Stefan mutual
diffusion matrix $D_{ij}$ (indices $1..m-1$, i.e. everything except the
polymer). The polymer block is handled separately — see
[`maxwell-stefan.md`](maxwell-stefan.md).

See [`nomenclature.md`](../nomenclature.md) for symbol definitions.

## 1. Infinite-dilution diffusion coefficient — Siddiqi & Lucas (1986)

Estimates $D_{ij}^{0}$, the diffusion coefficient of solute $i$ at infinite
dilution in solvent $j$, from the solvent's viscosity and both compounds'
molar volumes at the normal boiling point.

Aqueous solvent ($j$ = water):

$$
D_{ij}^{0} = 2.98\times10^{-7} \, \mu_{w}^{-1.026} \, V_{b,i}^{-0.5473} \, T
$$

Non-aqueous solvent:

$$
D_{ij}^{0} = 9.86\times10^{-8} \, \mu_{j}^{-0.907} \, V_{b,i}^{-0.45} \, V_{b,j}^{0.265} \, T
$$

- $\mu$ in mPa·s, $V_b$ in cm³/mol, $T$ in K → $D^0$ in cm²/s.
- Applies to non-polymer compound pairs only.

**Reference:** Siddiqi & Lucas (1986). See [`references.md`](../references.md).

## 2. Mutual diffusion from infinite-dilution values — Kooijman & Taylor (1991)

Combines the $D_{ij}^{0}$ matrix into composition-dependent mutual
diffusion coefficients using a logarithmic mixing rule:

$$
\ln D_{ij} = x_i \ln D_{ij}^{0} + x_j \ln D_{ji}^{0} + \sum_{\substack{k=1\\k\neq i,j}}^{m-1} \frac{x_k}{2}\left(\ln D_{ik}^{0} + \ln D_{jk}^{0}\right)
$$

for $i,j = 1, \dots, m-1$ (all non-polymer pairs).

**Reference:** Kooijman & Taylor (1991). See [`references.md`](../references.md).
