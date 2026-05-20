"""
Cobaya theory that provides the GSR (Generalized Slow-Roll)
primordial scalar power spectrum.

The basis tables (wkj, xkj) are loaded once in initialize().
Per-evaluation, mode amplitudes mj_0 .. mj_{n_basis-1} are
assembled from sampled parameters, combined with As/ns/n_run/n_runrun,
and the modified P_s(k) is computed on the basis k-grid, then
cubic-spline interpolated for arbitrary k queries.

YAML example
------------
theory:
  gsr_primordial_ps.GSRPrimordialPk:
    extra_args:
      file_wkj: /path/to/wkj.dat
      file_xkj: /path/to/xkj.dat
      k0: 0.05

params:
  As:
    prior: {min: 1.0e-9, max: 4.0e-9}
    ref: 2.1e-9
    latex: A_s
  ns:
    prior: {min: 0.9, max: 1.1}
    ref: 0.965
    latex: n_s
  mj_0:
    prior: {dist: norm, loc: 0, scale: 1}
    ref: 0.0
    latex: m_0
  mj_1:
    prior: {dist: norm, loc: 0, scale: 1}
    ref: 0.0
    latex: m_1
  # ... one mj_j per basis function in the wkj/xkj files
"""

import numpy as np
from scipy.interpolate import CubicSpline
from cobaya.theory import Theory
from cobaya.typing import empty_dict, InfoDict
from typing import Mapping


# ── basis I/O ──────────────────────────────────────────────────────────
def _read_basis_file(path: str):
  with open(path) as fh:
    header = fh.readline().split()
    n_k, n_basis = int(float(header[0])), int(float(header[1]))
    k = np.empty(n_k)
    vals = np.empty((n_k, n_basis))
    for i in range(n_k):
      row = fh.readline().split()
      k[i] = float(row[0])
      vals[i, :] = [float(x) for x in row[1 : n_basis + 1]]
  return k, vals


def _load_basis(file_wkj: str, file_xkj: str) -> dict:
  k_w, wkj = _read_basis_file(file_wkj)
  k_x, xkj = _read_basis_file(file_xkj)
  if wkj.shape != xkj.shape:
    raise ValueError(
      f"Incompatible shapes: wkj {wkj.shape} vs xkj {xkj.shape}"
    )
  if not np.allclose(k_w, k_x, rtol=1e-7):
    raise ValueError("k grids in wkj and xkj files are not consistent")
  return dict(k=k_w, wkj=wkj, xkj=xkj,
              n_k=len(k_w), n_basis=wkj.shape[1])


def _interp_at_pivot(k_grid, table, kpivot):
  n_basis = table.shape[1]
  vals = np.empty(n_basis)
  for j in range(n_basis):
    cs = CubicSpline(k_grid, table[:, j], bc_type="natural")
    vals[j] = cs(kpivot)
  return vals


# ── Cobaya theory ──────────────────────────────────────────────────────
class InPCReiPC(Theory):
  renames: Mapping[str, str] = empty_dict
  extra_args: InfoDict = {}

  def initialize(self):
    super().initialize()

    file_wkj = self.extra_args["file_wkj"]
    file_xkj = self.extra_args["file_xkj"]
    self._k0 = self.extra_args.get("k0", 0.05)

    basis = _load_basis(file_wkj, file_xkj)
    self._k = basis["k"]
    self._wkj = basis["wkj"]
    self._xkj = basis["xkj"]
    self._n_basis = basis["n_basis"]

    # precompute pivot splines (k0 is fixed across the chain)
    self._wj_pivot = _interp_at_pivot(self._k, self._wkj, self._k0)
    self._xj_pivot = _interp_at_pivot(self._k, self._xkj, self._k0)

    self.log.info(
      f"GSR basis loaded: n_k={len(self._k)}, "
      f"n_basis={self._n_basis}, k0={self._k0}"
    )

  def get_requirements(self):
    reqs = {"As": None, "ns": None}
    for j in range(self._n_basis):
      reqs[f"m{j}"] = None
    return reqs

  def get_can_provide(self):
    return ["primordial_scalar_pk"]

  def calculate(self, state, want_derived=True, **params_values_dict):
    As = self.provider.get_param("As")
    ns = self.provider.get_param("ns")

    # assemble mode amplitudes
    mj = np.array([
      self.provider.get_param(f"m{j}") for j in range(self._n_basis)
    ])

    # ── GSR computation (vectorised) ──────────────────────────────
    apivot = self._wj_pivot @ mj
    bpivot = self._xj_pivot @ mj

    I1SR2 = np.pi**2 * (1.0 - ns)**2 / 8.0

    I1CONTRIBUTION_PIVOT = (
      1.0 + 0.5 * (np.pi * (1.0 - ns) * bpivot + bpivot**2) / (1.0 + I1SR2)
    )

    a = self._wkj @ mj
    b = self._xkj @ mj

    I1CONTRIBUTION = (
      1.0 + 0.5 * (np.pi * (1.0 - ns) * b + b**2) / (1.0 + I1SR2)
    )

    lnrat = np.log(self._k / self._k0)
    ps_powerlaw = As * np.exp(lnrat * (ns - 1.0))

    ps = (
      ps_powerlaw
      * np.exp(a - apivot)
      * (I1CONTRIBUTION / I1CONTRIBUTION_PIVOT)
    )

    # build interpolator for arbitrary-k queries
    interp = CubicSpline(self._k, ps, bc_type="natural")

    state["primordial_scalar_pk"] = {
      "log_regular": False,
      "k": self._k,
      "Pk": ps,
      "effective_ns_for_nonlinear": ns,
    }

    arg_k = I1SR2 + 0.5 * (np.pi * (1.0 - ns) * b + b**2)
    I1PK = np.sqrt(max(np.max(arg_k), 0.0))

    state["I1PK"] = I1PK;
    state["derived"] = state.get("derived", {})
    state["derived"]["I1PK"] = I1PK

  def get_can_provide_params(self):
    return ['I1PK']

  def get_primordial_scalar_pk(self):
    return self.current_state["primordial_scalar_pk"]
