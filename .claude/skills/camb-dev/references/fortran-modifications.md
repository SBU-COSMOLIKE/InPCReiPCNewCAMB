# CAMB Fortran modifications: patterns and catalog

## The dispatch pattern (canonical example)

All model selection lives in `CAMB_ReadParams` in `camb.f90`, fenced and
driven by ini keys. This is the template for adding any new selectable
model class:

```fortran
    !VM: BEGINS - model dispatch: reionization
    block
      character(LEN=:), allocatable :: ReionModel, IPModel
      ReionModel = UpperCase(Ini%Read_String_Default('reionization_model', 'tanh'))
      if (ReionModel == 'BASIS') then
          if (allocated(P%Reion)) deallocate(P%Reion)
          allocate(TBasisReionization :: P%Reion)
      else if (ReionModel == 'EXP') then
          if (allocated(P%Reion)) deallocate(P%Reion)
          allocate(TExpReionization :: P%Reion)
      end if
      ! default (TANH) was set by CAMB_SetDefParams; no action needed
      call P%Reion%ReadParams(Ini)

      IPModel = UpperCase(Ini%Read_String_Default('initial_power_model', 'powerlaw'))
      if (IPModel == 'INFLATIONPC') then
          if (allocated(P%InitPower)) deallocate(P%InitPower)
          allocate(TInflationPCPower :: P%InitPower)
      end if
      ! default (POWERLAW) was set by CAMB_SetDefParams; no action needed
      call P%InitPower%ReadParams(Ini)
    end block
    !VM: ENDS
```

Corresponding `.ini` keys:

```ini
reionization_model = basis         ! or tanh (default) or exp
initial_power_model = inflationpc  ! or powerlaw (default)
```

Rules embedded in this pattern:

- The **default branch does nothing** — `CAMB_SetDefParams` already
  allocated the stock type, so the unmodified path is provably identical
  to upstream.
- `deallocate` before `allocate` of the polymorphic component; never
  assume the slot is empty.
- Dispatch is case-insensitive via `UpperCase`; ini values documented
  with their defaults in the comment and the key list.
- The required `use` additions (`use Reionization, only:
  TBasisReionization`, `use InitialPower, only: TInflationPCPower`) are
  part of the same fenced change set — a port is incomplete without them.

## Modification catalog (what exists and why)

### `reionization.f90` — `TBasisReionization`

Basis/principal-component x_e(z) histories as a new `TReionizationModel`
subtype. Stock `TTanhReionization` and `TExpReionization` untouched.
Pairs with the thermo and lsampling accuracy changes below — a basis
history with structure in x_e(z) is meaningless at default sampling.

### `InitialPower.f90` / `power_tilt_mj.f90` — `TInflationPCPower`

Mode-amplitude (mj) primordial power spectrum as a `TInitialPower`
subtype. Historically the Fortran route for GSR-type spectra
(`old_power_tilt_mj.f90` merged with `old_power_tilt.f90`); the
inference-pipeline route is now the Python `GSRPrimordialPk` Theory class
(see cobaya-integration.md), but the Fortran class is kept for ini-driven
standalone CAMB runs. Keep both consistent: same basis files, same
conventions, mj=0 recovers the same power law.

### `thermo.f90` / `ThermoData%inithermo` — nine changes

1. `nthermo = 30000` (denser thermodynamics grid).
2. Basis pre-pass block (evaluate the basis x_e before the main loop).
3. `fe` clamping (keep the free-electron fraction physical when basis
   amplitudes push it negative/above max).
4. Trapezoid `sdotmu` accumulation.
5. Trapezoid `actual_opt_depth`.
6. `tight_tau` fix (tight-coupling switch consistent with modified
   opacity history).
7. Opacity floor (numerical guard).
8. `our_derived_tau` assignment (derived optical depth stored for
   cross-checks and as a Cobaya derived parameter).
9. `SetRecombHistory` call placement.

These travel together. Porting a subset produces a build that runs and is
wrong; the trapezoid pair (4, 5) in particular must be consistent or the
derived tau disagrees with the integrated opacity.

### `lsampling.f90` — `initlval`, denser low-ell sampling

Six line changes, all inside the `if (CP%AccurateReionization)` branch:

| Old | New |
|---|---|
| `do lvar=11, 37,1` | `do lvar=11, 39,1` |
| `do lvar=11, 37,2` | `do lvar=11, 39,2` |
| `step = max(nint(5*Ascale),2)` | `step = 1` |
| `top=bot + step*10` | `top = 90` |
| `step=max(nint(20*Ascale),4)` | `step = 1` |
| `top=bot+step*2` | `top = 150` |
| `step=max(nint(25*Ascale),4)` | `step = 1` |
| `top=bot+step` | `top = 200` |

Net effect: every ell up to 200 sampled exactly (step = 1). Required for
reionization-basis EE features; without it the interpolated low-ell EE
smooths over exactly the structure the basis introduces.

### `CAMBparams` — added fields

```fortran
real(dl) :: our_derived_tau = 0._dl
```

Always with an explicit initializer. Lives wherever upstream currently
hosts the type (`model.f90` / `results.f90` — find it fresh each port,
upstream moves it).

## Porting discipline (lessons from the 2026 migration)

- Enumerate first: `grep -n "!VM" old_*.f90` before reading code. The
  block count is the port's definition of done.
- Old and new files are read **in full** before porting a block —
  upstream refactors relocate subroutines, and matching by surrounding
  context is the only safe method. Never apply a stored line-number diff.
- When the target file isn't in hand, emit a stand-alone patch document:
  full paste-ready blocks for structural changes, old → new tables for
  scattered one-liners (the lsampling table above is the format).
- A port session that hits a tool/length limit ends with a written
  "remaining work" list precise enough that the next session can resume
  without re-deriving anything.
- After porting, diff the unfenced remainder of each file against
  upstream — it must be empty. Accidental drift outside fences is the
  failure mode that makes the *next* port unbounded.

## Validation for Fortran changes

1. Default path: C_ell from the modified build with default ini keys is
   bit-identical to unmodified upstream CAMB.
2. Limit recovery: each custom model at its degenerate parameter point
   reproduces the stock model's C_ell (basis → tanh-equivalent; mj = 0 →
   power law).
3. Cross-build: each custom model agrees with the previous (old-CAMB)
   build to within the accuracy settings.
4. Derived quantities: `our_derived_tau` vs the input/intended tau as an
   internal consistency check on the opacity integration.
5. Accuracy: rerun with the accuracy knobs (nthermo, lAccuracyBoost,
   AccuracyBoost) increased one notch; the answer must be stable at the
   level the likelihoods care about.
