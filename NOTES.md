# Notes about the Free Volume Theory

## Generelized self-diffusion
The original selff-diffusion model of Vrentas&Duda is formulated for a binary system where the solvent penetrates into the polymer. In 1984, the same authors presented a a formualtion for the self diffussion coefficient $D_1$ and $D_2$ of the two penetrants in a polymer solvent system. 

$$
D_{2} = D_{02} \exp{\left(-\frac{w_{1}\hat{V}_{1}^{*}+w_{2}\hat{V}_{2}^{*}(\xi_{13}/\xi_{23}) + w_{3}\hat{V}_{3}^{*}\xi_{13}}{\hat{V}_{FH} / \gamma} \right)} 
$$

$$
D_{1} = D_{02} \exp{\left(-\frac{w_{1}\hat{V}_{1}^{*}(\xi_{23}/\xi_{13})+w_{2}\hat{V}_{2}^{*} + w_{3}\hat{V}_{3}^{*}\xi_{23}}{\hat{V}_{FH} / \gamma} \right)} 
$$

This expression were the base to Kubackza in its multicomponent expression: 

$$
D_{i} = D_{0i} \exp{\left({\frac{-E^{*}}{RT}}\right)} \exp{\left( -\frac{\sum_{j=1}^{m} w_{j}V_{j}^{*}\xi_{ip} / \xi_{jp}}{V_{FH} / \gamma} \right)} \quad i=1,2,...,m
$$



## About the $\xi$ parameter
The jumping unit ratio $\xi_{ip}$​ represents the ratio between the molar volume of the jumping unit of the diffusing species and that of the polymer. The molar jumping unit refers to the minimum volume required for a species to execute a molecular jump during diffusion.

$$
\xi_{ip} = \frac{\tilde{V}_{i}^{\circ}(0)}{\tilde{V}_{2j}^{*}}
$$

For solvents (or any small species diffusing within a polymer matrix), it is often reasonable to assume that the jumping unit volume corresponds to the molar volume at 0 K, which can be estimated using Bondi's group contribution method. On the other hand, estimating the polymer jumping unit volume is more complex. Hong (1995) proposed an empirical correlation for this parameter, which can be applied—with discretion—to well-characterized polymers such as polystyrene:

$$
\tilde{V}_{2j}^{*} (cm^3/mol) = 0.0925T_{g2} (K) + 69.47 \quad (T_{g2} < 295 K)
$$

$$
\tilde{V}_{2j}^{*} (cm^3/mol) = 0.6224T_{g2} (K) - 86.95 \quad (T_{g2} >= 295 K)
$$



## About the preexponetial factor $D_{0i}$
Kubaczka et al. (2014) employed the concept of friction coefficients from Bearman’s equation to develop a generalized expression for computing the mutual diffusion coefficient in polymer–solvent interactions. This formulation incorporates the self-diffusion coefficients of all components involved, including the polymer.

The pre-exponential factor of the solvent (or diffusing species) can be determined through non-linear regression of experimental density and viscosity data as a function of temperature. According to Hong et al. (1995), the pre-exponential factor D0D0​ is solvent-specific and independent of the polymer. In their analysis, they assume negligible energetic contributions (i.e., E=0E=0) and focus solely on the solvent's contribution to diffusion.

$$
\ln(\mu_{i})=\ln \left( \frac{0.124 \times 10^{-16} V_{ci}^{2/3} \rho_{i}RT}{M_{i}} \right)-\ln(D_{0i})+\frac{\hat{V}_{i}^{*}}{\frac{K_{1i}}{\gamma}(K_{2i}-T_{gi}+T)} \quad i = 1,2,3,... m -1
$$

In the foundational work of Vrentas and Duda (1977), the Free Volume Theory was introduced, providing expressions to compute the self-diffusion coefficients of both the solvent and the polymer. Later, Kubaczka et al. (2018) adapted these expressions to calculate the parameter DopDop​ using the Williams–Landel–Ferry (WLF) equation. However, this approach is limited to polymers for which WLF parameters are available.

The generalized form proposed by Kubaczka et al. (2018) is:

$$
\ln D_{0p} = \ln \left( \frac{\rho_{p}N_{a}}{36\eta_{p}} \right)\left( \frac{R^{2}}{M_{p}}\right) k_B T + \frac{K_{1p}}{\gamma}(K_{2p}+T-T_{gp})
$$

This expression captures the temperature dependence of the pre-exponential factor by incorporating viscosity, density, and free volume parameters of the polymer.

In this works we employ a polymer membrane where its principal component (POMS) has a lack of experimetnaal information avaible. So the parameter $D_{op}$ is fitting from experimental data. 

**Do somthing to incorporate a function that computes $D_{0p}$ for polymers with the WLF constatns**

## About the parameter $\hat{V}_{FH}$
Various authors (Hong, Zielinski, Vrentas & Duda) describe the general self-difussion in a binary solvent (1)-polymer(2) system like: 

$$
D_{1} = D_{01} \exp{\left({\frac{-E^{*}}{RT}}\right)} \exp{\left(\frac{-(w_{1}\hat{V}_{1}^{*}+w_{2}\xi\hat{V}_{2}^{*})}{\hat{V}_{FH} / \gamma} \right)} 
$$

Where the term in the denominator is defiend as: 
$$
\frac{\hat{V}_{FH}}{\gamma} =  w_{1} \frac{K_{11}}{\gamma} \left( K_{21} - T_{g1} + T  \right)  + w_{2}\frac{K_{12}}{\gamma}(K_{22}-T_{g2}+T)
$$

In 1997 Vrentas & Vrentas published a the methods to estimate the free volume parameters to compute the self-diffusion coefficient in solvet-polymer sistem. In this publication, they rewrite the definition of ${V_{FH}/}{\gamma}$ as:
$$
\frac{\hat{V}_{FH}}{\gamma} =  w_{1} \frac{K_{11}}{\gamma_1} \left( K_{21} - T_{g1} + T  \right)  + w_{2}\frac{\hat{V}_{FH2}}{\gamma_2}
$$

Now $\hat{V}_{FH2}$ depends on many other properties of the polymer. The authores recommend different expression for rubbery and glassy polymers. 

In 1988, Vrentas and Duda presented the influeuence of the Glass Transition on Solvent Self-Diffusion in Amorphous Polymers, here they explore the selfdiffusion in polymer solvent mixture below the glass transition temeprature. In this work apear the concepts like: specific volume of the equilibrum liquid polymer at the glass trnasition temeprature or the thermal expansion coefficient of the pure solvent. 

in 1994 they applied the free volume theory for solvent self.diffusion in polymer-solvent systems applied to rubbery mixture at temeprature above and below

*In this work*

For simplicity, this work employs the classic therminology to compute the denominator due the lack of information of the polymer that composes the membrane. For well known polymer with avaible data we recommend to employ the more complex approuch for better estimations. 















# References 
Vrentas, J. S., Duda, J. L., & Ling, H. ‐C. (1984). Self‐diffusion in polymer‐solvent‐solvent systems. Journal of Polymer Science: Polymer Physics Edition, 22(3), 459–469. https://doi.org/10.1002/pol.1984.180220308

Hong, S.-U. (1995). Prediction of Polymer/Solvent Diffusion Behavior Using Free-Volume Theory. Industrial & Engineering Chemistry Research, 34(7), 2536–2544. https://doi.org/10.1021/ie00046a040

Vrentas, J. S., & Duda, J. L. (1977). Diffusion in polymer–solvent systems. II. A predictive theory for the dependence of diffusion coefficients on temperature, concentration, and molecular weight. Journal of Polymer Science: Polymer Physics Edition, 15(3), 417–439. https://doi.org/10.1002/pol.1977.180150303

Kubaczka, A., Kamiński, W., & Marszałek, J. (2018). Predicting mass fluxes in the pervaporation process using Maxwell-Stefan diffusion coefficients. Journal of Membrane Science, 546, 111–119. https://doi.org/10.1016/j.memsci.2017.08.074

Vrentas, J. S., & Vrentas, C. M. (1998). Predictive methods for self-diffusion and mutual diffusion coefficients in polymer–solvent systems. European Polymer Journal, 34(5–6), 797–803. https://doi.org/10.1016/S0014-3057(97)00205-X
