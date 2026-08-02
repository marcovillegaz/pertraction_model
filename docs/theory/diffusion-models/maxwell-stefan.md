# Maxwell-Stefan formulation — polymer mutual diffusion, $[B]$, and $[\Gamma]$

Governs: `source/calculations/diffusion-models/kubaczka.m`,
`source/calculations/diffusion-models/maxwellStefan/Bmatrix.m`,
`source/calculations/diffusion-models/maxwellStefan/thermodynamicsFactors.m`.

Multicomponent mass transfer is modeled using the generalized
Maxwell-Stefan formulation, following the definitions in Taylor & Krishna
(1993) as adapted by Kubaczka for polymer-multicomponent systems. See
[`nomenclature.md`](../nomenclature.md) for symbol definitions. The
polymer is always the last component, index $m$.

## 1. Polymer mutual diffusion — Kubaczka (2014, 2018)

Non-polymer pairs $D_{ij}$ ($i,j < m$) come from
[`dilute-correlations.md`](dilute-correlations.md); polymer pairs
$D_{im}$ are derived here from the self-diffusion coefficients
$D_i^{self}$ (from [`free-volume-theory.md`](free-volume-theory.md)) using
Bearman's friction-coefficient formalism:

$$
K_{s,i} = 1 - \frac{x_i}{\displaystyle\sum_{j=1}^{m-1} x_j \dfrac{D_{i}^{self}}{D_{j}^{self}}}
$$

$$
\frac{1}{D_{i,m}} = \frac{1}{x_m}\left(\frac{K_{s,i}}{D_i^{self}} - \sum_{\substack{j=1\\j\neq i}}^{m-1} \frac{x_j}{D_{ij}}\right), \qquad D_{m,i} = D_{i,m}
$$

for $i = 1, \dots, m-1$.

## 2. Inverse drag matrix $[B]$

Assembles the full mutual diffusion matrix (non-polymer pairs from
§ dilute correlations, polymer pairs from § 1 above) into the
$(m-1)\times(m-1)$ inverse-drag matrix used in the Maxwell-Stefan flux
equations:

$$
B_{ii} = \frac{x_i}{D_{i,m}} + \sum_{\substack{k=1\\k\neq i}}^{m-1} \frac{x_k}{D_{ik}}
$$

$$
B_{ik} = -x_i\left(\frac{1}{D_{ik}} - \frac{1}{D_{i,m}}\right), \qquad k \neq i
$$

## 3. Thermodynamic factor matrix $[\Gamma]$ — Kubaczka (2014), eq. (17)

$$
\Gamma_{ik} = \delta_{ik}+x_{i}\left(\frac{\partial \ln \gamma_i}{\partial x_k} - \frac{\partial \ln \gamma_i}{\partial x_m} \right) \quad i,k = 1,2,\dots,m-1
$$

Because the activity-coefficient models (UNIFAC — see
[`../thermodynamic-models/unifac.md`](../thermodynamic-models/unifac.md))
are composed of many parameters and are difficult to differentiate
analytically, $\partial \ln\gamma_i/\partial x_k$ is approximated by
forward finite difference:

$$
\frac{\partial \ln\gamma_i}{\partial x_k} \approx \frac{\ln\gamma_i(x_k + h) - \ln\gamma_i(x_k)}{h}
$$

with step $h = 0.001$ in the current implementation
(`thermodynamicsFactors.m`). The activity-coefficient function returns a
vector of $\ln\gamma$ for all components, polymer last.

**References:** Kubaczka (2014); Kubaczka, Kamiński & Marszałek (2018);
Taylor & Krishna (1993). See [`references.md`](../references.md).
