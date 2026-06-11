# Cobaya–CAMB integration: the Python route

Preferred route for new physics that stock CAMB can ingest from outside:
a Cobaya `Theory` subclass computes the quantity in Python and hands it to
the unmodified `camb` theory block. Zero Fortran porting tax. The worked
example is `GSRPrimordialPk` (GSR primordial scalar P(k) from precomputed
basis tables), which replaced the Fortran `old_power_tilt_mj.f90` route in
the inference pipeline.

## The Theory-class pattern

- Subclass `cobaya.theory.Theory`, modeled on an existing block in the
  repo (the EuclidEmulator2 theory block is the in-house reference
  pattern).
- Load static inputs (basis tables `wkj`/`xkj`) once at initialization;
  precompute everything parameter-independent (pivot splines) there too.
  `calculate()` does only the parameter-dependent, vectorized math.
- Sampled parameters are the physics parameters (`As`, `ns`, `mj_0` ...
  `mj_{n-1}`). Parameters the model does not use are hardcoded to zero
  (`n_run = n_runrun = 0`), not exposed and ignored — the baseline is
  exactly `As * (k/k0)**(ns - 1)` with the GSR correction factor
  `exp(a - a_pivot) * (I1 / I1_pivot)` on top.
- Validate the degenerate limit numerically before wiring into the
  sampler: all `mj = 0` must recover the pure power law exactly.
- Input-file readers are defensive about the real files: skip trailing
  comment lines, assert the expected shape (420 k-points x 20 basis
  functions for the natural20 set), fail loudly on mismatch.

## The `external_primordial_pk` contract (hard-won)

To feed a custom primordial P(k) to CAMB through Cobaya:

- YAML: stock `camb:` theory block with
  `extra_args: {external_primordial_pk: True}`. Do **not** subclass the
  Cobaya CAMB wrapper for this — a `CAMBWithGSR` subclass was written and
  then deleted once it was established the flag is supported by the
  wrapper directly. Wrapper subclasses are a maintenance liability; use
  them only when the wrapper truly lacks a hook.
- The Theory class must publish a `primordial_scalar_pk` state dict with
  **exactly** the keys the wrapper consumes (found by reading
  `cobaya/theories/camb/camb.py` at the consumption site, ~line 684 in
  the version in use):

```python
state['primordial_scalar_pk'] = {
  'log_regular': False,
  'k': k,                              # 1D array
  'Pk': Pk,                            # same shape as k
  'effective_ns_for_nonlinear': ns_eff,
}
```

  With `log_regular: False` the wrapper calls `set_scalar_table(k, Pk)`.
  Missing keys fail with errors that point nowhere near the real cause.

The general lesson is the contract-discovery method: when a
Cobaya/CAMB integration fails, do not iterate on guesses. Open the wrapper
source at the exact consumption site, read which keys/arguments it
requires, and conform. Two guessed fixes failed; one source read fixed it.

## Established Theory modules in the pipeline

- `GSRPrimordialPk` — GSR primordial scalar spectrum (above).
- `InPCReiPC` — initial-power + reionization PC theory module used with
  the Planck 2018 likelihood stack (plik_lite TTTEEE, lowl TT, lowl EE
  sroll2) and the custom I1C likelihood in the profile-likelihood /
  emcee+MPI minimizer.

New modules follow the same conventions so they compose in one YAML.

## Practical integration notes

- Cobaya patch maintenance: in-house patches to cobaya must track
  upstream's formatting (black reformatting, clik → clipy migration have
  both broken patches before). Keep patches minimal and re-anchor them
  after upstream formatting churn.
- Suppress known-noisy likelihood logging through Cobaya's logging tree,
  not Python warnings:
  `logging.getLogger("planck_2018_lowl.ee_sroll2").setLevel(logging.ERROR)`
  (`warnings.filterwarnings` does not catch Cobaya's logging output).
- When a sampled parameter keeps rejecting at the prior box (the
  `theta_MC_100` case: reference value ~0.7 sigma from the prior lower
  bound), the symptom is a bimodal evaluation-time distribution across
  MPI ranks (fast prior rejections vs full CAMB evaluations). Check the
  prior geometry before suspecting the code.
- Standalone scripts in this layer follow the group Python style:
  2-space indentation, concise production scripts, explicit key-ordered
  parameter handling (never rely on dict `.values()` ordering for a
  parameter vector).
