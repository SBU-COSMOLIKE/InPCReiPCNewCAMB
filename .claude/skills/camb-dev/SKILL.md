---
name: camb-dev
description: Methodology for modifying CAMB (Fortran Boltzmann code) and integrating it with Cobaya inside the Cocoa ecosystem, Miranda group style. Use this skill whenever working with CAMB in any way — writing or porting Fortran modifications (reionization, initial power, thermodynamics, dark energy models), migrating !VM-fenced changes to a new CAMB release, building Cobaya Theory classes that feed CAMB (primordial P(k), background quantities), implementing background dark-sector models like the axion phantom-mirage scenario, or debugging CAMB/Cobaya integration errors. Also trigger on mentions of reionization basis, TBasisReionization, TInflationPCPower, external_primordial_pk, thermo.f90, lsampling, w(a) tables, or effective-fluid dark energy in a CAMB context.
---

# CAMB Development (Miranda group methodology)

This skill encodes how custom physics is added to CAMB and kept alive across
upstream releases in the Cocoa ecosystem: VM-fenced Fortran modifications
(`TBasisReionization`, `TInflationPCPower`, thermo/lsampling changes),
Python-side Cobaya `Theory` classes, and background dark-sector models
(axion phantom-mirage). The style is conservative by design — CAMB feeds
Planck-grade likelihoods, and a silent physics regression is worse than any
amount of porting tedium.

Read `references/fortran-modifications.md` before touching any `.f90` file.

Read `references/cobaya-integration.md` before writing or debugging a Cobaya
`Theory` class or CAMB `extra_args`.

## Core principles

1. **Fence everything.** Every modification to upstream CAMB Fortran is
   wrapped in `!VM BEGINS` / `!VM ENDS` markers (`!VM: BEGINS` variant also
   in use — grep for `!VM`). The fences are the contract that makes
   migration to new CAMB releases mechanical: `grep -n "!VM"` enumerates
   the complete modification surface. Code outside fences is upstream and
   must stay byte-identical to upstream.
2. **New physics is a new class behind an ini key; defaults are untouched.**
   Add polymorphic types (`TBasisReionization` alongside stock
   `TTanhReionization`/`TExpReionization`; `TInflationPCPower` alongside
   `TInitialPowerLaw`) and select them via dispatch on an ini string
   (`reionization_model`, `initial_power_model`). The stock default path
   must produce bit-identical output to unmodified CAMB. Never edit the
   behavior of an existing model class to add physics.
3. **Limit recovery is the primary validation.** Every new model must
   exactly reduce to a known model at a parameter limit, and that limit is
   tested numerically, not assumed: mode amplitudes mj=0 recovers the pure
   power law; basis reionization with the tanh-equivalent coefficients
   recovers tanh C_ell; axion fraction f_ax → 0 recovers LambdaCDM. A new
   model without a tested limit is not done.
4. **Read the integration layer's source, don't guess its API.** Cobaya's
   CAMB wrapper (`camb.py`) defines undocumented contracts (e.g. the
   `primordial_scalar_pk` state dict keys). When an integration fails,
   open the wrapper source at the call site and read what it actually
   expects. This resolved the GSR integration in one step after two
   guessed attempts failed.
5. **Python before Fortran.** If the new physics is a quantity CAMB can
   ingest from outside — a primordial P(k) table, a w(a)/background table,
   a derived prior — implement it as a Cobaya `Theory` class in Python and
   feed stock CAMB. Reserve Fortran for physics that must live inside the
   Boltzmann/thermodynamics evolution (reionization history shapes,
   recombination, opacity). Fortran mods cost a porting tax at every CAMB
   release; Python classes don't.
6. **Independent cross-checks in Mathematica.** Background quantities —
   growth integrals D(a) ∝ H(a) ∫ da'/(a' H(a'))^3, comoving distances,
   w(a) solutions, Limber mappings — are verified independently via the
   Wolfram MCP before trusting a pipeline number. The check is a different
   implementation, not a re-run of the same code.
7. **Accuracy settings travel with the physics.** Models with sharp
   features (basis reionization, localized w(a) transitions) need the
   matching sampling-density changes (denser ell sampling in
   `lsampling.f90`, `nthermo = 30000`, `AccurateReionization`). Shipping
   the model without its accuracy settings produces silently biased C_ell.

## Migration workflow (new CAMB release)

This is the procedure used for the 2026 CAMB port; follow it exactly.

1. Collect old (modified) and new (upstream) files side by side; prefix old
   ones `old_` so both can coexist in one directory.
2. `grep -n "!VM"` every old file to enumerate the modification surface
   before reading anything else. Count blocks; the port is done when every
   block is accounted for in the new tree.
3. Read each old file and its new counterpart fully — upstream refactors
   move code, so never port by line number, always by surrounding context.
4. Port block by block. Where upstream changed the surrounding structure,
   adapt the block to the new structure but keep the physics identical and
   the fence markers intact.
5. For files not yet available, produce a stand-alone patch document: the
   exact code block to paste, plus an old → new table for scattered
   one-line changes (the `lsampling.f90` six-line table is the model).
6. Validate by limit recovery and by comparing C_ell of the default path
   against unmodified upstream CAMB (must be bit-identical) and of each
   custom model against the old build (must agree to the old build's
   accuracy settings).

## Style

- Fortran inside fences matches the **surrounding upstream CAMB style**
  (indentation, `_dl` kinds, CamelCase type names with `T` prefix,
  `block`/`end block` for local scopes in dispatch code). The fence carries
  ownership; style mismatches just make future ports harder.
- Patches against upstream Python (cobaya, CAMB's python layer) match
  upstream formatting exactly — black-formatted where upstream is black —
  to keep maintained diffs minimal across upstream reformatting.
- Standalone group-owned Python (Theory classes, analysis scripts) uses
  2-space indentation and the group's concise, production-script style: no
  speculative generality, no format options nobody asked for.
- New CAMBparams fields get explicit initializers
  (`real(dl) :: our_derived_tau = 0._dl`).
- Hardcode to zero what the model doesn't use (e.g. `n_run`, `n_runrun` in
  the GSR baseline) instead of carrying dead parameters through the math.
- Ini keys are lowercase snake_case strings dispatched case-insensitively
  (`UpperCase(Ini%Read_String_Default('reionization_model', 'tanh'))`).

## Map of the modified tree

- `camb.f90` — `CAMB_ReadParams`: model-selection dispatch block
  (deallocate/allocate polymorphic `P%Reion`, `P%InitPower`, then
  `ReadParams(Ini)`); `use` statements for the custom modules.
- `reionization.f90` — `TBasisReionization` (PC/basis x_e(z) histories).
- `InitialPower.f90` / `power_tilt_mj.f90` — `TInflationPCPower`
  (mode-amplitude primordial spectrum; superseded for GSR by the Python
  route, kept for ini-driven runs).
- `thermo.f90` (`ThermoData` / `inithermo`) — `nthermo = 30000`, basis
  pre-pass, `fe` clamping, trapezoid `sdotmu` / `actual_opt_depth`,
  `tight_tau` fix, opacity floor, `our_derived_tau` assignment.
- `lsampling.f90` (`initlval`) — denser ell sampling under
  `AccurateReionization` (see the patch table in the Fortran reference).
- `CAMBparams` (model/results module) — added derived/storage fields
  (e.g. `our_derived_tau`).
- Python side — Cobaya `Theory` classes (e.g. `GSRPrimordialPk`,
  `InPCReiPC`) feeding stock `camb:` via `extra_args`.

## Review checklist for CAMB patches

- [ ] All changes inside `!VM` fences; `grep "!VM"` block count matches
      the patch description.
- [ ] Default model path bit-identical to upstream (C_ell diff is zero).
- [ ] Limit-recovery test exists, runs, and passes for the new model.
- [ ] No behavior change to stock classes; new physics in new types only.
- [ ] Accuracy settings (ell sampling, nthermo) shipped with the model
      that needs them.
- [ ] Integration contracts (state dict keys, extra_args) verified against
      the wrapper source, with the source line cited in the PR.
- [ ] Background numbers cross-checked independently (Mathematica) when
      the patch introduces or alters H(z), w(a), or growth.
- [ ] Style matches surroundings; ini keys documented with defaults.
