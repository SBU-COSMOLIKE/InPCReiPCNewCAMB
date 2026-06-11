# Axion / phantom-mirage background models in the CAMB pipeline

Context for implementing and forecasting background dark-sector models —
the science driver behind the (beta, p) program in the DOE/NSF EC and
Roman proposals.

## The model class

DESI DR2 favors dynamical dark energy at ~3 sigma, with an expansion that
looks like Lambda at low and high redshift and dynamical only at
intermediate times. A single-fluid w(z) reproducing this requires phantom
crossing, which healthy particle-physics models tend to forbid. The
group's alternative class ("phantom mirage"): the apparent phantom
behavior is a mis-specification artifact. Concrete realizations under
study:

1. Multi-field dark energy (Hu 2005).
2. **Lambda plus an axion** with mass tuned so the field is frozen at
   early times (w ~ -1, dark-energy-like) and oscillates today
   (cycle-averaged w ~ 0, CDM-like), with the transition near z_c ~ 0.5
   (Liu:2025bss, Caldwell:2025inn).
3. Explicit dark-sector couplings (Miranda:2017rdk, Khoury:2025txd).

Data constrain only the total expansion history, not the partition of the
dark sector — so a w0waCDM fit to (2) infers phantom crossing that isn't
there. The falsification observable is growth: SPCA (smooth paradigm of
cosmic acceleration) predicts G(z) fully determined by H(z) and Omega_m;
phantom-mirage models break that link.

## The (beta, p) detection statistic

```
G(a | w0, wa, Omega_m, beta, p) / G(a | w0, wa, Omega_m)
    = 1 - beta * ( Omega_DE(a) / (1 - Omega_m) )^p
```

beta = 0 is the sharp SPCA null. Anchors (verify before reuse — these are
proposal-era values): Lin et al. 2023 prefer beta = 0.149 +/- 0.052 under
HaloFit-based nonlinear modeling; preliminary Roman cosmic-shear Fisher
gives >= 3 sigma detection if |beta| ~ 0.1. Framing discipline from the
referee rounds: (beta, p) is a *detector* of SPCA departure, not a
parameterization of any microphysics; a nonzero beta falsifies SPCA, a
null is consistent but not a validation.

Known limitation to carry into any implementation: the axion transition
at z_c ~ 0.5 produces a *localized* feature in G(z) that the smooth
power law may absorb only partially. The agreed robustness check is to
repeat with a more flexible G(z) basis (principal components, in the
spirit of Mortonson, Hu & Huterer 2009) and confirm the beta inference is
stable.

## Implementation routes (in order of preference)

1. **Python, background-only (default).** Solve the model background
   externally; feed CAMB an effective dark-energy w(a) table through the
   Python interface (`DarkEnergyFluid` + `set_w_a_table(a, w)` via Cobaya
   `extra_args` / a Theory class). The effective single-fluid w(a)
   combines Lambda + axion: w_eff = (p_Lambda + p_ax)/(rho_Lambda +
   rho_ax). No Fortran, no porting tax. This is the route consistent with
   camb-dev principle 5.
2. **Mathematica for the model side.** The background solve (Klein-Gordon
   for the field, or the cycle-averaged effective fluid once H < m), the
   growth ODE, and all cross-check integrals
   (D(a) ∝ H(a) ∫ da'/(a' H(a'))^3, chi(z), Limber k = (l+1/2)/chi) run
   in Mathematica via the Wolfram MCP. This is both the prototyping tool
   and the independent verification of whatever lands in the pipeline.
3. **Fortran `TDarkEnergyModel` subtype** only if perturbation-level
   treatment inside CAMB becomes necessary (it has not so far: the model
   class is chosen precisely because it leaves the dark matter
   perturbative equations unchanged, which also keeps the EFTofLSS
   kernels valid for scale-independent growth changes — a claim flagged
   in review as needing explicit verification, so verify it for any new
   variant).

The (beta, p) growth modification itself is applied downstream
(CosmoLike / EFTofLSS side) as a scale-independent rescaling of G(a); it
does not require a CAMB modification.

## Model → (beta, p) mapping pipeline

1. Solve the model background for H(a) (axion fraction f_ax, mass m or
   z_c, plus Lambda).
2. Solve growth in the true model: the axion sources matter clustering
   after z_c (CDM-like) and not before — the source term changes across
   the transition.
3. Fit w0waCDM to the model's expansion history (the mirage fit) and
   compute the SPCA reference G(a | w0, wa, Omega_m).
4. Form the ratio and fit the (beta, p) template; record the best-fit
   locus and the template residual (the residual is the quantitative
   version of the "localized feature" concern above).
5. Place the locus against the Fisher/forecast reach and the Lin et al.
   anchor. Never quote a model's (beta, p) without also quoting the
   template residual.

Every step gets an independent Mathematica check before its output is
trusted in a proposal figure or pipeline input.

## June 2026 axion session — TO BE FILLED BY VIVIAN

The week of June 8, 2026 included a working session on the axion model
(alongside the NSF EC port) whose concrete outputs are not reproduced
here: the solver code, the implementation route actually chosen, and any
computed (beta, p) loci or sigma(beta) staging numbers. Do not invent
these values. If this reference is being used and those results matter,
ask Vivian for the session's notebook/outputs and append them to this
section verbatim: model parameters (f_ax, m or z_c), the computed
(beta, p) per parameter point, template residuals, and forecast numbers.
