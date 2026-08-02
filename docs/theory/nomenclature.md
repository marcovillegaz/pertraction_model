# Nomenclature

Single source of truth for symbols used across `docs/theory/`. Every theory
doc should use these symbols and link back here instead of redefining them
locally. If a doc needs a symbol not listed here, add it here first.

Convention: index `1..m-1` are the mobile species (solvent/solute), index
`m` is always the polymer (membrane) when a system includes one.

## Composition

| Symbol | Meaning | Units |
|---|---|---|
| $m$ | number of components in the system (last index = polymer, when present) | – |
| $x_i$ | molar fraction of component $i$ | – |
| $w_i$ | mass fraction of component $i$ | – |
| $n$ | degree of polymerization (UNIFAC-FV input) | – |

## Diffusion coefficients

| Symbol | Meaning | Units |
|---|---|---|
| $D_{ij}^{0}$ | mutual diffusion coefficient at infinite dilution ($i$ dilute in $j$) | cm²/s |
| $D_{ij}$ | mutual (Maxwell-Stefan) diffusion coefficient | cm²/s |
| $D_i^{self}$ | self-diffusion coefficient of component $i$ | cm²/s |
| $D_{0i}$ | pre-exponential (reference) diffusion coefficient of $i$ in the FVT expression | cm²/s |
| $[B]$ | inverse drag coefficient matrix, $(m-1)\times(m-1)$ | s/cm² |
| $[\text{Đ}]$ | effective Fick diffusivity matrix, $[B]^{-1}[\Gamma]$ | cm²/s |
| $\mu$ | dynamic viscosity | mPa·s (correlation), Pa·s (property store) |
| $V_b$ | molar volume at the normal boiling point | cm³/mol |
| $T$ | absolute temperature | K |

## Free volume theory (Vrentas & Vrentas / Kubaczka)

| Symbol | Meaning | Units |
|---|---|---|
| $\xi_{ip}$ | jumping-unit ratio of species $i$ relative to reference $p$, $\xi_{ip}=V_i^{*}/V_p^{*}$ | – |
| $V_i^{*}$ | specific critical hole free volume of $i$ required for a jump | cm³/g |
| $\hat{V}_{FH}/\gamma$ | hole free volume of the mixture, weighted average over components | cm³/g |
| $K_{1i}/\gamma$ | free-volume parameter ($\beta_i$ in code) | cm³/(g·K) |
| $K_{2i}$ | free-volume parameter (temperature offset) | K |
| $T_{gi}$ | glass transition temperature of $i$ | K |
| $E^{*}$ | activation energy for a diffusive jump | cal/mol |
| $R$ | universal gas constant (1.987 cal/mol·K in this codebase) | cal/(mol·K) |

## Thermodynamics (UNIFAC / Maxwell-Stefan)

| Symbol | Meaning | Units |
|---|---|---|
| $\gamma_i$ | activity coefficient of component $i$ | – |
| $[\Gamma]$ | thermodynamic factor matrix, $(m-1)\times(m-1)$ | – |
| $\delta_{ik}$ | Kronecker delta | – |
| $r_i, q_i$ | UNIFAC volume and surface-area parameters of compound $i$ | – |
| $R_k, Q_k$ | UNIFAC group volume and surface-area parameters | – |
| $\theta_i, \phi_i$ | UNIFAC surface-area and volume fractions | – |
| $a_{mn}$ | UNIFAC group-group interaction parameter | K |
| $\psi_{mn}$ | UNIFAC interaction term, $\exp(-a_{mn}/T)$ | – |
| $z$ | UNIFAC coordination number (= 10) | – |

## Where symbols are defined precisely

Some symbols carry a slightly different definition per model family (e.g.
$D_{ij}$ means something different in the dilute correlations vs. the
polymer mutual-diffusion step). When that happens, the nuance is explained
in the relevant model doc — this table gives the general meaning and units
so notation stays consistent across documents.
