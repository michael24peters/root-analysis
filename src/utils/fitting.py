import math
import numpy as np
from scipy.special import erf
from scipy.integrate import cumulative_trapezoid


def _trapz(y, x):
    """
    Simple trapezoidal integration: ∫ y dx
    """
    y = np.asarray(y)
    x = np.asarray(x)
    dx = np.diff(x)
    avg = 0.5 * (y[:-1] + y[1:])
    return np.sum(avg * dx)

# ─── Histogram loader ─────────────────────────────────────────────────────────

def load_histogram(root_path, branch="tag_dtf_m", tree="tree",
                   xmin=480.0, xmax=620.0, nbins=80):
    """
    Read a ROOT TTree branch with uproot and return a binned histogram.
    Flattens jagged arrays (multiple candidates per event) before binning,
    using the same logic as src/matplotlib/plot_mass.py:30.

    Returns
    -------
    centers : ndarray (nbins,)   bin centre positions [MeV]
    counts  : ndarray (nbins,)   observed counts per bin (float)
    errors  : ndarray (nbins,)   per-bin uncertainty: max(sqrt(counts), 1)
    """
    import uproot
    with uproot.open(root_path) as f:
        arr = f[tree][branch].array(library="np")
    flat = np.concatenate(arr) if arr.dtype == object else arr.ravel()
    counts, edges = np.histogram(flat, bins=nbins, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    errors  = np.maximum(np.sqrt(counts.astype(float)), 1.0)
    return centers, counts.astype(float), errors


# ─── Chebyshev background ─────────────────────────────────────────────────────
"""
The Chebyshev polynomials are defined by:
T0(x̃) = 1
T1(x̃) = x̃
T2(x̃) = 2x̃² - 1

where

x̃ = (2x - (xmin + xmax)) / (xmax - xmin)  maps x ∈ [xmin, xmax] to x̃ ∈ [-1, 1].
Parameterizing this way makes the variable normalized, centered at zero,
and symmetric, which means T0 and T2 are even functions and T1 is odd. This
means <T0, T1> = <T1, T2> = 0, so the coefficients c0 and c1 are uncorrelated in
a fit. <T0, T2> ≠ 0, so c1 will be correlated with the overall normalization 
(c0=0) but not with the slope (c0≠0).

We use a 2nd-order Chebyshev polynomial for the background shape, which is
B(x̃) = 1 + c0·T1(x̃) + c1·T2(x̃)
"""

def _xtilde(x, xmin, xmax):
    """
    Define parameterized, normalized [-1, 1] variable for the Chebyshev
    polynomials.

    x̃ = (2x - (xmin + xmax)) / (xmax - xmin)
    """
    return (2.0 * np.asarray(x, float) - (xmin + xmax)) / (xmax - xmin)


def _cheb_norm(xmin, xmax, c1):
    """
    Analytic normalization ∫_{xmin}^{xmax} B(x̃) dx over [xmin, xmax].
    T1(x̃) = x̃ is odd, so its integral over [-1, 1] vanishes.
    T2(x̃) = 2x̃² - 1 contributes: (xmax-xmin)/2 · (2 - 4c1/3) = (xmax-xmin)·(1 - c1/3).
    """
    return (xmax - xmin) * (1.0 - c1 / 3.0)


def cheb_bkg_pdf(x, xmin, xmax, c0, c1):
    """
    Normalized 2nd-order Chebyshev background PDF [1/MeV].
    B(x̃) = 1 + c0·T1(x̃) + c1·T2(x̃), where T1(x̃) = x̃, T2(x̃) = 2x̃² - 1.
    Matches RooFit's RooChebychev convention.
    """
    u = _xtilde(x, xmin, xmax)
    return (1.0 + c0 * u + c1 * (2.0 * u**2 - 1.0)) / _cheb_norm(xmin, xmax, c1)


def _cheb_antideriv(u, c0, c1):
    """Indefinite integral of B(x̃) = 1 + c0·x̃ + c1·(2x̃² - 1) w.r.t. x̃."""
    return u + c0 * u**2 / 2.0 + c1 * (2.0 * u**3 / 3.0 - u)


def cheb_cdf(x, xmin, xmax, c0, c1):
    """
    Analytic Chebyshev CDF normalized to [0, 1] over [xmin, xmax].
    Computed via antiderivative; equals 0 at xmin and 1 at xmax.
    """
    norm  = _cheb_norm(xmin, xmax, c1)
    scale = (xmax - xmin) / 2.0    # Jacobian: dx = (xmax-xmin)/2 · du
    u_lo  = -1.0                    # xtilde(xmin) = -1
    u     = _xtilde(x, xmin, xmax)
    return scale * (_cheb_antideriv(u, c0, c1)
                    - _cheb_antideriv(u_lo, c0, c1)) / norm


# ─── Gaussian signal ──────────────────────────────────────────────────────────
"""
f(x) = 1/(σ√(2π)) · exp(-0.5 * ((x - μ)/σ)²)
"""

def gauss_pdf(x, mean, sigma, xmin, xmax):
    """
    Gaussian signal PDF [1/MeV], normalized over [xmin, xmax].
    Uses scipy.special.erf for finite-range normalization.

    $$
    \int_a^b G(x)dx=\frac{1}{\sigma\sqrt{2\pi}}\int_a^b e^{-\frac{(x-\mu)^2}{2\sigma^2}}dx
    =\frac{1}{2}\left[\text{erf}\left(\frac{b-\mu}{\sigma\sqrt{2}}\right)-\text{erf}\left(\frac{a-\mu}{\sigma\sqrt{2}}\right)\right]
    $$

    ∫{a,b} ​G(x)dx = 1/2 ​[erf((b - mu) / (sigma * sqrt2)) − erf((a - mu) / (sigma * sqrt2))]
    """
    x      = np.asarray(x, float)
    sqrt2  = math.sqrt(2.0)
    norm   = 0.5 * (erf((xmax - mean) / (sigma * sqrt2))
                    - erf((xmin - mean) / (sigma * sqrt2)))
    raw    = np.exp(-0.5 * ((x - mean) / sigma)**2) / (sigma * math.sqrt(2.0 * math.pi))
    return raw / norm


def gauss_cdf(x, mean, sigma, xmin, xmax):
    """
    Gaussian CDF normalized to [0, 1] over [xmin, xmax].
    """
    x      = np.asarray(x, float)
    sqrt2  = math.sqrt(2.0)
    erf_lo = erf((xmin - mean) / (sigma * sqrt2))
    erf_hi = erf((xmax - mean) / (sigma * sqrt2))
    erf_x  = erf((x    - mean) / (sigma * sqrt2))
    norm   = 0.5 * (erf_hi - erf_lo)
    return 0.5 * (erf_x - erf_lo) / norm


# ─── Double Crystal Ball signal ───────────────────────────────────────────────
"""
A Gaussian core with independent power-law tails on each side of the peak:

          ╔  A_L · (B_L − t)^{−n_L}    if   t < −α_L       (left  tail)
  f(t) =  ║  exp(−t²/2)                if  −α_L ≤ t ≤ α_R  (Gaussian core)
          ╚  A_R · (B_R + t)^{−n_R}    if   t >  α_R        (right tail)

  where   t = (x − μ) / σ
  and the continuity + differentiability constants are
          A_{L,R} = (n / |α|)^n  · exp(−α²/2)
          B_{L,R} =  n / |α|    −  |α|

α_{L,R} > 0 : distance from the mean (in σ units) where the PDF switches
              from Gaussian to power-law.
n_{L,R} > 1 : power-law exponent; larger → harder / more abrupt tail falloff.
"""

def _dcb_tail_constants(alpha, n):
    """
    Continuity constants A, B for one DCB tail (alpha > 0).
    At the transition point t = -|alpha| (left) or t = |alpha| (right):
      A = (n/|α|)^n · exp(-α²/2)
      B = n/|α| - |α|
    """
    abs_a = abs(alpha)
    A = (n / abs_a) ** n * math.exp(-0.5 * alpha * alpha)
    B = n / abs_a - abs_a
    return A, B


def dcb_unnorm(x, mean, sigma, alphaL, nL, alphaR, nR):
    """
    Unnormalized Double Crystal Ball PDF.
    """
    x  = np.asarray(x, float)
    t  = (x - mean) / sigma
    AL, BL = _dcb_tail_constants(alphaL, nL)
    AR, BR = _dcb_tail_constants(alphaR, nR)
    out = np.empty_like(t)
    mc = (-alphaL <= t) & (t <= alphaR)
    ml = t < -alphaL
    mr = t >  alphaR
    out[mc] = np.exp(-0.5 * t[mc]**2)
    out[ml] = AL * (BL - t[ml]) ** (-nL)
    out[mr] = AR * (BR + t[mr]) ** (-nR)
    return out


def _dcb_grid(xmin, xmax, mean, sigma, alphaL, nL, alphaR, nR, n_grid=1000):
    """
    Evaluate unnormalized DCB on a dense grid for CDF/norm computation.
    """
    xg = np.linspace(xmin, xmax, n_grid)
    fg = dcb_unnorm(xg, mean, sigma, alphaL, nL, alphaR, nR)
    return xg, fg


def dcb_pdf(x, mean, sigma, alphaL, nL, alphaR, nR, xmin, xmax):
    """
    Normalized DCB PDF over [xmin, xmax].

    The DCB has no analytical integrla, so it is instead computed numerically
    via trapezoid rule (integral as sum of trapezoids) on a dense grid, i.e., 
    many slices.
    """
    # Evaluate DCB on evenly spaced grid
    xg, fg = _dcb_grid(xmin, xmax, mean, sigma, alphaL, nL, alphaR, nR)
    # Integrate via trapezoid rule
    norm   = _trapz(fg, xg)
    return dcb_unnorm(np.asarray(x, float),
                      mean, sigma, alphaL, nL, alphaR, nR) / norm


def dcb_cdf(x, mean, sigma, alphaL, nL, alphaR, nR, xmin, xmax):
    """
    DCB CDF normalized to [0, 1] over [xmin, xmax].
    Computed via cumulative trapezoid on a dense grid, then interpolated.
    """
    # Evaluate DCB on evenly spaced grid
    xg, fg   = _dcb_grid(xmin, xmax, mean, sigma, alphaL, nL, alphaR, nR)
    cdf_grid = np.zeros_like(fg)
    # This implementation is O(N). Using _trapz in a loop would be O(N^2).
    # Width of each trapezoid: dx_i = xg[i+1] - xg[i]
    # Avg height of each trapezoid: avg_i = 0.5 * (fg[i] + fg[i+1])
    # Area of each trapezoid: avg_i * dx_i
    # Running total: cdf_grid[i+1] = cdf_grid[i] + avg_i * dx_i
    cdf_grid[1:] = np.cumsum(0.5 * (fg[:-1] + fg[1:]) * np.diff(xg))
    # Normalize to [0, 1]
    cdf_grid = cdf_grid / cdf_grid[-1]
    return np.interp(np.asarray(x, float), xg, cdf_grid)


# ─── iminuit cost function factories ──────────────────────────────────────────
#
# Each factory returns an iminuit.cost.ExtendedBinnedNLL instance.
# The inner `scaled_cdf(xe, *params)` function must return:
#   (n_total, cdf)  where n_total = n_sig + n_bkg  and
#   cdf[i] = expected events with x < xe[i], so cdf[0]=0, cdf[-1]=n_total.
# iminuit infers parameter names from the inner function's signature.
#

def make_gauss_cost(counts, bin_edges, xmin, xmax):
    """
    ExtendedBinnedNLL cost for Gaussian signal + 2nd-order Chebyshev background.
    Free parameters: mean, sigma, n_sig, n_bkg, c0, c1
    """
    from iminuit.cost import ExtendedBinnedNLL

    def scaled_cdf(xe, mean, sigma, n_sig, n_bkg, c0, c1):
        return (
            gauss_cdf(xe, mean, sigma, xmin, xmax) * n_sig
            + cheb_cdf(xe, xmin, xmax, c0, c1) * n_bkg
        )

    return ExtendedBinnedNLL(counts, bin_edges, scaled_cdf)


def make_dcb_cost(counts, bin_edges, xmin, xmax):
    """
    ExtendedBinnedNLL cost for DCB signal + 2nd-order Chebyshev background.
    Free parameters: mean, sigma, alphaL, nL, alphaR, nR, n_sig, n_bkg, c0, c1
    """
    from iminuit.cost import ExtendedBinnedNLL

    def scaled_cdf(xe, mean, sigma, alphaL, nL, alphaR, nR, n_sig, n_bkg, c0, c1):
        return (
            dcb_cdf(xe, mean, sigma, alphaL, nL, alphaR, nR, xmin, xmax) * n_sig
            + cheb_cdf(xe, xmin, xmax, c0, c1) * n_bkg
        )

    return ExtendedBinnedNLL(counts, bin_edges, scaled_cdf)


def make_dcb_sym_cost(counts, bin_edges, xmin, xmax):
    """
    ExtendedBinnedNLL cost for symmetric DCB (α_R ≡ α_L, n_R ≡ n_L) + Chebyshev.
    Free parameters: mean, sigma, alpha, n, n_sig, n_bkg, c0, c1  (8 vs 10)
    """
    from iminuit.cost import ExtendedBinnedNLL

    def scaled_cdf(xe, mean, sigma, alpha, n, n_sig, n_bkg, c0, c1):
        return (
            dcb_cdf(xe, mean, sigma, alpha, n, alpha, n, xmin, xmax) * n_sig
            + cheb_cdf(xe, xmin, xmax, c0, c1) * n_bkg
        )

    return ExtendedBinnedNLL(counts, bin_edges, scaled_cdf)


# ─── iminuit fit wrapper ───────────────────────────────────────────────────────

def run_fit(cost, initial_values, limits):
    """
    Run MIGRAD (minimisation) then HESSE (covariance / symmetric errors).

    Parameters
    ----------
    cost           : iminuit cost function (e.g. ExtendedBinnedNLL)
    initial_values : dict {param_name: start_value}
    limits         : dict {param_name: (lo, hi)}

    Returns
    -------
    m : iminuit.Minuit
        Fitted Minuit object. Key attributes:
          m.values     – best-fit parameter values (MnUserParameters-like dict)
          m.errors     – HESSE symmetric errors
          m.valid      – True if MIGRAD converged
          m.fval       – minimum function value (2 * NLL)
          m.nfit       – number of free floating parameters
          m.fmin.edm   – estimated distance to minimum
    """
    from iminuit import Minuit
    m = Minuit(cost, **initial_values)
    for param, lim in limits.items():
        m.limits[param] = lim
    m.migrad()
    if m.valid:
        m.hesse()
    return m
