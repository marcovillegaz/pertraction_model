# UNIFAC model for computation of activity coefficients

Governs: `source/calculations/thermodynamic-models/unifac/*.m` (shared
terms), `unifacLegacy.m`, `unifacFV.m`, `unifacVdwFV.m` (variant
composition), and `source/classes/ThermoModel/*.m` (class layer + factory).

## Activity coefficient

In a mixture of ordinary liquids, it is customary to define an activity
coefficient $\gamma_i$ for component $i$:

$$
\gamma_{i} = \frac{a_i}{x_i}
$$

where $a$ is the activity and $x$ the mole fraction. The activity is
defined as the ratio of the fugacity of component $i$ in the mixture to
that in a standard state, usually chosen to be pure liquid $i$ at the
system temperature and pressure.

UNIFAC is a group-contribution method for estimating $\gamma_i$ from group
interaction parameters rather than pure-component vapor-liquid data. Three
selectable variants are implemented, all built on the same combinatorial
($C$) and residual ($R$) terms and differing only in whether/how a
free-volume ($FV$) correction is added:

$$
\ln\gamma_i = \ln\gamma_i^{C} + \ln\gamma_i^{R} + \ln\gamma_i^{FV}
$$

| Variant | Terms | Class | Function | `RequiredCompoundProps` |
|---|---|---|---|---|
| Legacy UNIFAC | $C + R$ | `UNIFACLegacyModel` | `unifacLegacy.m` | *(none)* |
| UNIFAC-FV | $C + R + FV$ | `UNIFACFVModel` | `unifacFV.m` | `FVP`, `density` |
| UNIFAC-vdW-FV | $C + R + FV_{vdW}$ | `UNIFACvdWFVModel` | `unifacVdwFV.m` | `FVP`, `density` |

Select a variant via the factory rather than constructing a class directly:

```matlab
model = ThermoModel.create("unifac-fv", compoundsLib, unifacLib);   % or "unifac", "unifac-vdw-fv"
LnGamma = model.computeActivityCoefficient(T, x);
```

`ThermoModel.create` calls `validateCompoundData()` before returning, which
checks every compound in the library against the chosen variant's
`RequiredCompoundProps` and raises `ThermoModel:MissingCompoundData` if a
required field is absent *or present but empty* (e.g. an `FVP` row with no
value) -- see [`../nomenclature.md`](../nomenclature.md) for what each FVP
row means.

Consumed by
[`../diffusion-models/maxwell-stefan.md`](../diffusion-models/maxwell-stefan.md)
§ 3 to build $[\Gamma]$, via `thermodynamicsFactors(thermoModel, T, x)`,
which takes a `ThermoModel` object directly.

---

## Legacy UNIFAC

Fredenslund, Jones & Prausnitz (1975). See [`../references.md`](../references.md).

$\ln\gamma_i = \ln\gamma_i^{C} + \ln\gamma_i^{R}$ -- no free-volume term.
Suitable for non-polymer systems. `UNIFACLegacyModel` /
`unifacLegacy.m`.

### Combinatorial part -- `unifacCombinatorial.m`

$$
r_k = \sum_{\text{groups}} R_{\text{group}}\,\nu_{k,\text{group}}, \qquad q_k = \sum_{\text{groups}} Q_{\text{group}}\,\nu_{k,\text{group}}
$$

$$
l_k = \frac{z}{2}(r_k - q_k) - (r_k - 1), \qquad z = 10
$$

For a polymer component, $r$ is scaled by the degree of polymerization $n$:
$r_{\text{polymer}} \leftarrow n \cdot r_{\text{polymer}}$ (done once, in
`unifacGroupParams.m`, shared by every term that reads $r$). Legacy UNIFAC
calls `unifacGroupParams` with $n=1$ -- a no-op, matching its lack of a
polymer correction.

$$
\phi_k = \frac{r_k x_k}{\sum_l r_l x_l}, \qquad \theta_k = \frac{q_k x_k}{\sum_l q_l x_l}
$$

$$
\ln\gamma_k^{C} = \ln\frac{\phi_k}{x_k} + \frac{z}{2}q_k\ln\frac{\theta_k}{\phi_k} - \frac{\phi_k}{x_k}\sum_l x_l l_l + l_k
$$

### Residual part (group contributions) -- `unifacResidual.m`

Identical for all three variants -- reused by UNIFAC-FV and UNIFAC-vdW-FV
below rather than duplicated.

Group interaction term: $\psi_{mn} = \exp(-a_{mn}/T)$.

Residual activity coefficient of group $k$ in a reference solution of
pure compound $i$ ($X_{k,i}$ = mole fraction of group $k$ in pure $i$,
$\theta_{m,i}$ = area fraction of group $m$ in pure $i$):

$$
\ln\Gamma_k^{(i)} = Q_k\left[1 - \ln\left(\sum_m \theta_{m,i}\psi_{mk}\right) - \sum_m \frac{\theta_{m,i}\psi_{km}}{\sum_n \theta_{n,i}\psi_{nm}}\right]
$$

The same expression, with $X_k$ and $\theta_m$ computed over the actual
mixture instead of the pure compound, gives $\ln\Gamma_k$ (mixture).

$$
\ln\gamma_i^{R} = \sum_{k} \nu_{k,i}\left(\ln\Gamma_k - \ln\Gamma_k^{(i)}\right)
$$

where $\nu_{k,i}$ is the number of groups $k$ in compound $i$.

---

## UNIFAC-FV

Oishi & Prausnitz (1978). See [`../references.md`](../references.md).

$\ln\gamma_i = \ln\gamma_i^{C} + \ln\gamma_i^{R} + \ln\gamma_i^{FV}$, reusing
the combinatorial and residual terms above and adding a free-volume
correction for polymer-solvent systems. `UNIFACFVModel` / `unifacFV.m`.

### Free-volume part -- `unifacFreeVolume.m`

Hard-core volume $v_h = 15.17\,r$ (cm³/mol, Bondi-type); free volume
$v_{fv} = v - v_h$, where $v_k = MW_k/\rho_k(T)$ is the molar volume from
the compound's density function ($r$ here already includes the polymer
scaling from the combinatorial term above).

$$
\phi_k^{fv} = \frac{x_k v_{fv,k}}{\sum_l x_l v_{fv,l}}, \qquad \phi_k^{h} = \frac{x_k v_{h,k}}{\sum_l x_l v_{h,l}}
$$

$$
\ln\gamma_k^{FV} = \ln\frac{\phi_k^{fv}}{\phi_k^{h}} + \frac{\phi_k^{h}-\phi_k^{fv}}{x_k}
$$

---

## UNIFAC-vdW-FV

Kannan, Duda & Danner (2005), "A free-volume term based on the van der
Waals partition function for the UNIFAC model". See
[`../references.md`](../references.md).

$\ln\gamma_i = \ln\gamma_i^{C} + \ln\gamma_i^{R} + \ln\gamma_i^{FV_{vdW}}$,
reusing the combinatorial and residual terms above with a free-volume term
based on the van der Waals partition function, distinct from UNIFAC-FV's
Bondi hard-core-volume term. `UNIFACvdWFVModel` / `unifacVdwFV.m`.

### Free-volume part (van der Waals) -- **NEEDS FORMULA**

The citation above is now known, but the van der Waals-partition-function
free-volume expression itself has not yet been extracted from Kannan,
Duda & Danner (2005) into this document.

**Current implementation state:** `unifacVdwFV.m` reuses
`unifacFreeVolume.m` (the Bondi hard-core-volume term from UNIFAC-FV)
unchanged as a placeholder, so `UNIFACvdWFVModel` is currently numerically
**identical** to `UNIFACFVModel`. It is safe to call (tested, produces
finite output) but is **not yet a distinct, validated model** -- do not
treat its output as representing the van der Waals correction.

**To resolve:** extract the free-volume equation(s) from Kannan, Duda &
Danner (2005) into this section, then implement them as a new function
(e.g. `unifacFreeVolumeVdw.m`) alongside `unifacFreeVolume.m`, and wire it
into `unifacVdwFV.m` in place of the TODO-marked call.

---

## Notes

- Compounds with $x_i = 0$ have their $\ln\gamma_i$ set to 0 (no defined
  activity coefficient at zero concentration).
- The degree of polymerization $n$ is a constructor option on
  `UNIFACFVModel`/`UNIFACvdWFVModel` (default 1000, matching `main.m`'s
  prior hardcoded value) or `ThermoModel.create(..., "PolymerizationDegree", n)`.
  Legacy UNIFAC does not use it.

**References:** Fredenslund, Jones & Prausnitz (1975); Oishi & Prausnitz
(1978); Kannan, Duda & Danner (2005); Poling, Prausnitz & O'Connell (2001),
ch. 8.10. See [`../references.md`](../references.md).
