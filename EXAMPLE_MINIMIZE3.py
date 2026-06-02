import warnings, os, psutil, sys
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings(
    "ignore",
    message=".*column is deprecated.*",
    module=r"sacc\.sacc"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*invalid value encountered.*"
)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r".*overflow encountered*"
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=r".*Hartlap correction*"
)
import copy, argparse, random, time, emcee 
import numpy as np
from cobaya.yaml import yaml_load
from cobaya.model import get_model
from getdist import IniFile
from schwimmbad import MPIPool
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
_affinity_set = False

def enforce_affinity():
    local_rank = int(os.environ.get("OMPI_COMM_WORLD_LOCAL_RANK",
                     os.environ.get("MPI_LOCALRANKID", 0)))
    omp_threads = int(os.environ.get("OMP_NUM_THREADS", 1))
    first_core = local_rank * omp_threads
    last_core  = first_core + omp_threads - 1
    try:
        p = psutil.Process()
        p.cpu_affinity(list(range(first_core, last_core + 1)))
    except Exception as e:
        rank = int(os.environ.get("OMPI_COMM_WORLD_RANK", 0))
        print(f"[Rank {rank}] Failed to set affinity: {e}")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
parser = argparse.ArgumentParser(prog='EXAMPLE_MINIMIZE_LCDM_CAMB1')
parser.add_argument("--nstw",
                    dest="nstw",
                    help="Number of likelihood evaluations (steps) per temperature per walker",
                    type=int,
                    default=200)
parser.add_argument("--root",
                    dest="root",
                    help="Name of the Output File",
                    default="./projects/lsst_y1/")
parser.add_argument("--outroot",
                    dest="outroot",
                    help="Name of the Output File",
                    default="example_min1")
# need to use parse_known_args because of mpifuture 
args, unknown = parser.parse_known_args() 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
yaml_string = r"""
likelihood:
  planck_2018_highl_plik.TTTEEE_lite: 
    path: ./external_modules/
    clik_file: plc_3.0/hi_l/plik_lite/plik_lite_v22_TTTEEE.clik
  planck_2018_lowl.TT: 
    path: ./external_modules
  planck_2018_lowl.EE_sroll2: null 
  I1C:
    external: "lambda _self: -np.inf if _self.provider.get_param('I1PK') > 0.7071067811 else 0.0"
    requires:
      I1PK: null
params:
  logA:
    prior:
      min: 2.8
      max: 3.2
    ref:
      dist: norm
      loc: 3.0448
      scale: 0.05
    proposal: 0.05
    latex: \log(10^{10} A_\mathrm{s})
    drop: true
  As:
    value: 'lambda logA: 1e-10*np.exp(logA)'
    latex: A_\mathrm{s}
  ns: 
    prior:
      min: 0.92
      max: 1.1
    ref:
      dist: norm
      loc: 0.96605
      scale: 0.005
    proposal: 0.005
    latex: n_\mathrm{s}
  theta_MC_100:
    prior:
      min: 0.99
      max: 1.5
    ref:
      dist: norm
      loc: 1.04109
      scale: 0.0015
    proposal: 0.0015
    latex: 100\theta_\mathrm{MC}
    drop: true
    renames: theta
  cosmomc_theta:
    value: 'lambda theta_MC_100: 1.e-2*theta_MC_100'
    derived: false
  H0:
    latex: H_0
    min: 20
    max: 100
  omegabh2:
    prior:
      min: 0.01
      max: 0.03
    ref:
      dist: norm
      loc: 0.022383
      scale: 0.005
    proposal: 0.005
    latex: \Omega_\mathrm{b} h^2
  omegach2:
    prior:
      min: 0.08
      max: 0.16
    ref:
      dist: norm
      loc: 0.12011
      scale: 0.01
    proposal: 0.01
    latex: \Omega_\mathrm{c} h^2
  tau:
    prior:
      min: 0.03
      max: 0.11
    ref:
      dist: norm
      loc: 0.055
      scale: 0.01
    proposal: 0.005
    latex: \tau_\mathrm{reio}
  w:
    value: -1
    latex: w_{0,\mathrm{DE}}
  w0pwa:
    value: -1
    latex: w_{0,\mathrm{DE}}+w_{a,\mathrm{DE}}
    drop: true
  wa:
    value: 'lambda w0pwa, w: w0pwa - w'
    latex: w_{a,\mathrm{DE}}
  mnu:
    value: 0.06
  omegam:
    derived: true
  age:
    derived: true
    latex: '{\rm{Age}}/\mathrm{Gyr}'
  rdrag:
    derived: true
    latex: r_\mathrm{drag}
  yheused:
    derived: true
    latex: Y_P^\mathrm{BBN}
  omegan2:
    derived: true
    latex: \Omega_\mathrm{\\nu} h^2
  zstar:
    derived: true
    latex: z_*
  rstar:
    latex: r_*
  thetastar:
    latex: 100\theta_*
  DAstar:
    latex: D_\mathrm{A}/\mathrm{Gpc}
  zdrag:
    latex: z_\mathrm{drag}
  kd:
    latex: k_\mathrm{D}
  thetad:
    latex: 100\theta_\mathrm{D}
  zeq:
    latex: z_\mathrm{eq}
  m0:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_0
  m1:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_1
  m2:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_2
  m3:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_3
  m4:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_4
  m5:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_5
  m6:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_6
  m7:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_7
  m8:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{8}
  m9:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{9}
  m10:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{10}
  m11:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{11}
  m12:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{12}
  m13:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{13}
  m14:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{14}
  m15:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{15}
  m16:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{16}
  m17:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{17}
  m18:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{18}
  m19:
    prior:
      min: -5 
      max: 5
    ref:
      dist: norm
      loc: 0.0
      scale: 0.05
    proposal: 0.05
    latex: m_{19}
theory:
  InPCReiPC:
    provides: ['I1PK']
    extra_args:
      file_wkj: ./projects/InPCReiPC/data/wj_natural20smin200smax2e4.dat
      file_xkj: ./projects/InPCReiPC/data/xj_natural20smin200smax2e4.dat
      k0: 0.05
  camb:
    path: ./external_modules/code/CAMB
    stop_at_error: False
    use_renames: True
    external_primordial_pk: True
    extra_args:
      halofit_version: takahashi
      dark_energy_model: ppf
      AccuracyBoost: 1.0  
"""
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
model = get_model(yaml_load(yaml_string))
import logging
logging.getLogger("planck_2018_lowl.ee_sroll2").setLevel(logging.ERROR)
def chi2(p):
    #p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    if isinstance(p, dict):
      p = [float(p[k]) for k in model.parameterization.sampled_params()]
    if np.any(np.isinf(p)) or  np.any(np.isnan(p)):
      raise ValueError(f"At least one parameter value was infinite (CoCoa) param = {p}")
    point = dict(zip(model.parameterization.sampled_params(), p))
    res1 = model.logprior(point,
                          make_finite=False)
    if np.isinf(res1) or np.any(np.isnan(res1)):
      return 1.e20
    res2 = model.loglike(point,
                         make_finite=False,
                         cached=False,
                         return_derived=False)
    if np.isinf(res2) or np.isnan(res2):
      return 1.e20
    return -2.0*(res1+res2)
def chi2v2(p):
    p = [float(v) for v in p.values()] if isinstance(p, dict) else p
    point = dict(zip(model.parameterization.sampled_params(), p))
    logposterior = model.logposterior(point, as_dict=True)
    chi2likes=-2*np.array(list(logposterior["loglikes"].values()))
    chi2prior=-2*np.atleast_1d(model.logprior(point,make_finite=False))
    return np.concatenate((chi2likes, chi2prior))
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def min_chi2(x0, 
             cov, 
             fixed=-1, 
             nstw=200,
             nwalkers=5,
             pool=None):
    def mychi2(params, *args):
        z, fixed, T = args
        params = np.array(params, dtype='float64')
        if fixed > -1:
            params = np.insert(params, fixed, z)
        return chi2(p=params)/T

    if fixed > -1:
        z      = x0[fixed]
        x0     = np.delete(x0, (fixed))
        args = (z, fixed, 1.0)
        
        cov = np.delete(cov, (fixed), axis=0)
        cov = np.delete(cov, (fixed), axis=1)
    else:
        args = (0.0, -2.0, 1.0)

    def logprob(params, *args):
        global _affinity_set
        if not _affinity_set:
          _affinity_set = True
          start_time = time.time()
          res = mychi2(params, *args)
          etime = time.time() - start_time
          rank = int(os.environ.get("OMPI_COMM_WORLD_RANK",0))
          print(f"Emcee: Like Eval Time: {etime:.4f} secs and MPI Rank: {rank}")
        else:
          res = mychi2(params, *args)
        if (res > 1.e19 or np.isinf(res) or  np.isnan(res)):
          return -np.inf
        else:
          return -0.5*res
    
    class GaussianStep:
       def __init__(self, stepsize=0.2):
           self.cov = stepsize*cov
       def __call__(self, x):
           return np.random.multivariate_normal(x, self.cov, size=1)   
    
    ndim        = int(x0.shape[0])
    nwalkers    = int(nwalkers)
    nstw        = int(nstw)
    temperature = np.array([1.0, 0.25, 0.1, 0.025, 0.005, 0.001], dtype='float64')
    ntemp       = len(temperature)
    stepsz      = temperature/4.0

    start_time = time.time()
    mychi2(GaussianStep(stepsize=0.001)(x0)[0,:], *args)
    elapsed_time = time.time() - start_time

    # emcee's DEMove/DESnookerMove use the red-blue split — within one step, 
    # the red half is updated using positions from the blue half (parallel), 
    # then the blue half is updated using the new red positions (parallel). 
    # Peak parallelism is nwalkers/2, not nwalkers.
    pool_size = pool.size if hasattr(pool, 'size') else pool.comm.Get_size() - 1
    batches_per_step = 2 * int(np.ceil(nwalkers / 2 / max(pool_size, 1)))
    wall_hours = elapsed_time * nstw * ntemp * batches_per_step / 3600.
    print(f"Estimated wall time: {wall_hours:.4f} hours "
          f"(nwalkers={nwalkers}, pool_size={pool_size}, "
          f"batches/step={batches_per_step})")

    partial_samples = []
    partial = []
    for i in range(len(temperature)):
        x = [] # Initial point
        for j in range(nwalkers):
            x.append(GaussianStep(stepsize=stepsz[i])(x0)[0,:])  
        sampler = emcee.EnsembleSampler(nwalkers, 
                                        ndim, 
                                        logprob, 
                                        args=(args[0], args[1], temperature[i]),
                                        moves=[(emcee.moves.DEMove(), 0.8),
                                               (emcee.moves.DESnookerMove(), 0.2)],
                                        pool=pool)    
        sampler.run_mcmc(np.array(x, dtype='float64'), 
                         nstw, 
                         skip_initial_state_check=True)
        samples = sampler.get_chain(flat=True, discard=0)
        j = np.argmin(-1.0*np.array(sampler.get_log_prob(flat=True)))
        partial_samples.append(samples[j])
        partial.append(mychi2(samples[j], *args))
        x0 = copy.deepcopy(samples[j])
        sampler.reset()    
        j = np.argmin(np.array(partial))
        print(f"Partial ({i+1}/{len(temperature)}): "
              f"params = {partial_samples[j]}, and chi2 = {partial[j]}")
    # min chi2 from the entire emcee runs
    j = np.argmin(np.array(partial))
    result = [partial_samples[j], partial[j]]
    return result
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
def prf(x0, nstw, cov, fixed=-1, nwalkers=5, pool=None):
    res =  min_chi2(x0=np.array(x0, dtype='float64'), 
                    fixed=fixed,
                    cov=cov, 
                    nstw=nstw, 
                    nwalkers=nwalkers,
                    pool=pool)
    return res
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    enforce_affinity() # enforce affinity (so Hybrid MPI-OpenMP works)!
    with MPIPool() as pool:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        dim = model.prior.d()  
        n_workers = max(pool.comm.Get_size() - 1, 1) # master doesn't compute
        nwalkers  = max(3*dim, 2*n_workers)
        nwalkers += nwalkers % 2                     # force even for red/blue split
        nstw = args.nstw
        (x0, results) = model.get_valid_point(max_tries=50, 
                                              ignore_fixed_ref=False,
                                              logposterior_as_dict=True)
        # 1st: Get covariance --------------------------------------------------
        cov = model.prior.covmat(ignore_external=True) # cov from prior
        
        # 2nd: Run Procoli -----------------------------------------------------
        res = np.array(list(prf(np.array(x0, dtype='float64'), 
                               fixed=-1, 
                               nstw=nstw,
                               nwalkers=nwalkers,
                               pool=pool,
                               cov=cov)), dtype="object")
        xf = np.array([res[0]],dtype='float64')
        
        # 3rd Append derived parameters ----------------------------------------
        xf = np.column_stack((xf, 
                              np.array([chi2v2(d) for d in xf], dtype='float64'),
                              res[1]))
        # 4th Save output file -------------------------------------------------
        os.makedirs(os.path.dirname(f"{args.root}chains/"), exist_ok=True)
        names = list(model.parameterization.sampled_params().keys()) # Cobaya Call
        hd = names 
        hd = hd + list(model.info()['likelihood'].keys()) + ["prior"] + ["chi2"]
        np.savetxt(f"{args.root}chains/{args.outroot}.txt", 
                   xf,
                   fmt="%.12e",
                   header=f"nswt (evals/Temp/walker)={nstw}\n"+' '.join(hd),
                   comments="# ")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------