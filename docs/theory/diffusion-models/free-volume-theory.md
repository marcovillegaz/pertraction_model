# Free Volume Theory — self-diffusion in polymer-solvent systems

Governs: `source/calculations/diffusion-models/vrentasVrentas.m`.

See [`nomenclature.md`](../nomenclature.md) for symbol definitions.

## Background

The original Vrentas & Duda self-diffusion model is formulated for a
binary solvent-polymer system. Vrentas, Duda & Ling (1984) extended it to
$D_1$ and $D_2$, the self-diffusion coefficients of two penetrants (1, 2)
in a polymer (3):

$$
D_{2} = D_{02} \exp{\left(-\frac{w_{1}\hat{V}_{1}^{*}+w_{2}\hat{V}_{2}^{*}(\xi_{13}/\xi_{23}) + w_{3}\hat{V}_{3}^{*}\xi_{13}}{\hat{V}_{FH} / \gamma} \right)}
$$

$$
D_{1} = D_{02} \exp{\left(-\frac{w_{1}\hat{V}_{1}^{*}(\xi_{23}/\xi_{13})+w_{2}\hat{V}_{2}^{*} + w_{3}\hat{V}_{3}^{*}\xi_{23}}{\hat{V}_{FH} / \gamma} \right)}
$$

Kubaczka generalized this to $m$ mobile components in a multicomponent
polymer system — this is the form implemented in `vrentasVrentas.m`:

$$
D_{i} = D_{0i} \exp{\left({\frac{-E^{*}}{RT}}\right)} \exp{\left( -\frac{\sum_{j=1}^{m} w_{j}V_{j}^{*}\,\xi_{ip}/\xi_{jp}}{\hat{V}_{FH} / \gamma} \right)} \quad i = 1,2,\dots,m
$$

> **Implementation check:** `vrentasVrentas.m` currently computes the
> activation term as `exp(-E(i)/R*T)`. In MATLAB `/` and `*` share
> precedence and evaluate left-to-right, so this evaluates as
> $-\frac{E_i}{R}\cdot T$, not $-\frac{E_i}{RT}$ as in the formula above.
> Verify against this doc when the API-mismatch fix (migration plan,
> Priority 2a) is applied.

## The jumping-unit ratio $\xi_{ip}$

$\xi_{ip}$ is the ratio between the molar volume of the jumping unit of
the diffusing species and that of the reference species $p$ (the polymer,
in this codebase):

$$
\xi_{ip} = \frac{\tilde{V}_{i}^{\circ}(0)}{\tilde{V}_{p}^{*}}
$$

For small diffusing species, the jumping-unit volume is commonly taken as
the molar volume at 0 K (estimable via Bondi's group contribution method).
Estimating the polymer's jumping-unit volume is harder. Hong (1995)
proposed an empirical correlation, applicable with discretion to
well-characterized polymers such as polystyrene:

$$
\tilde{V}_{p}^{*} \; (\text{cm}^3/\text{mol}) = 0.0925\,T_{g,p}\,(\text{K}) + 69.47 \quad (T_{g,p} < 295\,\text{K})
$$

$$
\tilde{V}_{p}^{*} \; (\text{cm}^3/\text{mol}) = 0.6224\,T_{g,p}\,(\text{K}) - 86.95 \quad (T_{g,p} \geq 295\,\text{K})
$$

In `vrentasVrentas.m`, $\xi_i$ is simplified to $V_i^{*}/V_m^{*}$ where $m$
is the last (polymer) index — i.e. it is computed directly from the
per-compound `FVP` hole-volume column rather than from the Hong
correlation above.

## The pre-exponential factor $D_{0i}$

Kubaczka et al. (2014) used Bearman's friction-coefficient formalism to
generalize the mutual diffusion coefficient for polymer-solvent
interactions, incorporating the self-diffusion coefficients of all
components (including the polymer).

For solvents, $D_{0i}$ can be regressed from experimental density and
viscosity vs. temperature, assuming negligible energetic contribution
($E \approx 0$, Hong et al. 1995):

$$
\ln(\mu_{i})=\ln \left( \frac{0.124 \times 10^{-16} V_{ci}^{2/3} \rho_{i}RT}{M_{i}} \right)-\ln(D_{0i})+\frac{\hat{V}_{i}^{*}}{\frac{K_{1i}}{\gamma}(K_{2i}-T_{gi}+T)} \quad i = 1,\dots,m-1
$$

For the polymer, Kubaczka et al. (2018) adapted this via the
Williams–Landel–Ferry (WLF) equation:

$$
\ln D_{0p} = \ln \left( \frac{\rho_{p}N_{a}}{36\eta_{p}} \right)\left( \frac{R^{2}}{M_{p}}\right) k_B T + \frac{K_{1p}}{\gamma}(K_{2p}+T-T_{gp})
$$

This is limited to polymers with known WLF parameters. **In this work**,
the membrane polymer (POMS) lacks sufficient experimental data, so
$D_{0p}$ is instead obtained by fitting to experimental data rather than
computed from WLF constants.

> **Open item:** add a function that computes $D_{0p}$ for polymers with
> known WLF constants, for use when a well-characterized polymer replaces
> POMS.

## The hole free-volume parameter $\hat{V}_{FH}/\gamma$

The general binary solvent (1)–polymer (2) self-diffusion expression
(Hong, Zielinski, Vrentas & Duda):

$$
D_{1} = D_{01} \exp{\left({\frac{-E^{*}}{RT}}\right)} \exp{\left(\frac{-(w_{1}\hat{V}_{1}^{*}+w_{2}\xi\hat{V}_{2}^{*})}{\hat{V}_{FH} / \gamma} \right)}
$$

with the denominator classically defined as:

$$
\frac{\hat{V}_{FH}}{\gamma} =  w_{1} \frac{K_{11}}{\gamma} \left( K_{21} - T_{g1} + T  \right)  + w_{2}\frac{K_{12}}{\gamma}(K_{22}-T_{g2}+T)
$$

Vrentas & Vrentas (1998) rewrote this in terms of the polymer's own
hole free volume:

$$
\frac{\hat{V}_{FH}}{\gamma} =  w_{1} \frac{K_{11}}{\gamma_1} \left( K_{21} - T_{g1} + T  \right)  + w_{2}\frac{\hat{V}_{FH2}}{\gamma_2}
$$

$\hat{V}_{FH2}$ depends on further polymer properties, with different
recommended expressions for rubbery vs. glassy polymers (see Vrentas &
Duda 1988 for the sub-$T_g$ treatment, and their 1994 extension to rubbery
mixtures above and below $T_g$).

**In this work**, given the lack of polymer (POMS) data, the classic
(1977) denominator form is used directly — the generalized multicomponent
version implemented in `vrentasVrentas.m`:

$$
\frac{\hat{V}_{FH}}{\gamma} = \sum_{j=1}^{m} w_{j}\,\beta_{j}\,(\delta_{j} + T), \qquad \beta_j = \frac{K_{1j}}{\gamma}, \quad \delta_j = K_{2j} - T_{gj}
$$

For a polymer with well-characterized free-volume data, prefer the more
complex Vrentas & Vrentas (1997/1998) approach instead.

**References:** Vrentas & Duda (1977); Vrentas, Duda & Ling (1984); Hong
(1995); Vrentas & Vrentas (1998); Kubaczka et al. (2018). See
[`references.md`](../references.md).
