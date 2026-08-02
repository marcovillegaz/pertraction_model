# Pure-component property models

Governs: `source/calculations/property-models/GCVOL60.m`,
`source/calculations/property-models/sastriRao.m`,
`source/calculations/property-models/TynCalus.m`.

Group-contribution correlations used by `loadCompoundData.m` to build the
per-compound density and viscosity function handles stored in
`CompoundsLibrary`, and by the diffusion correlations for molar volume at
the boiling point. See [`../nomenclature.md`](../nomenclature.md) for
symbol definitions.

## Density — GCVOL (Ihmels, 2003)

Group-contribution estimate of liquid density as a function of
temperature, using Elbro group parameters $(A, B, C, n)$ per group:

$$
\rho(T) = \frac{1000 \cdot MW}{\displaystyle\sum_{\text{groups}} n\,(A + BT + CT^2)}
$$

$\rho$ in kg/m³, $T$ in K, $MW$ in g/mol.

## Viscosity — Sastri & Rao (1992)

Group-contribution estimate of liquid viscosity as a function of
temperature, using group parameters $(\Delta\mu, \Delta N)$ and normal
boiling point $T_b$:

$$
\mu_b = \sum_{\text{groups}} \Delta\mu \cdot \text{occurrence}, \qquad N = 0.2 + \sum_{\text{groups}} \Delta N \cdot \text{occurrence}
$$

$$
\mu(T) = 0.001\,\mu_b \left[\exp\left((4.5396 + 1.0309\ln T_b)\left(1 - \frac{(3-2T_r)^{0.19}}{T_r} - 0.38\ln(T_r)(3-2T_r)^{-0.81}\right)\right)\right]^{-N}
$$

where $T_r = T/T_b$. Output in Pa·s.

## Molar volume at boiling point — Tyn & Calus (1975)

$$
V_b = 0.285\,V_c^{1.048}
$$

$V_b, V_c$ in cm³/mol. Used as the $V_b$ input to the Siddiqi & Lucas
correlation — see
[`../diffusion-models/dilute-correlations.md`](../diffusion-models/dilute-correlations.md).

**References:** Ihmels (2003); Sastri & Rao (1992); Tyn & Calus (1975).
See [`../references.md`](../references.md).
