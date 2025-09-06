import math

def q_asimov(f, L, a, b, c):
    """
    Asimov test statistic for a single Poisson bin (no backgrounds):
      q(f) = 2 [ mu(f) - mu0 + mu0 * ln(mu0 / mu(f)) ],
    with mu(f) = L * (c + a f + b f^2), mu0 = L * c.
    Units: f in TeV^-4, a,b,c in fb, L in fb^-1.
    """
    mu0 = L * c
    muf = L * (c + a*f + b*f*f)
    if muf <= 0:
        return float('inf')  # outside physical region
    return 2.0 * (muf - mu0 + mu0 * math.log(mu0 / muf))

def _limit_one_side(L, a, b, c, q_target=3.84, side="+"):
    """
    Find f on the chosen side of 0 such that q(f)=q_target.
    Robust outward geometric scan + bisection.
    """
    # start from f=0 and march outward
    step = 1.0 if side == "+" else -1.0
    f_prev = 0.0
    f = step
    for _ in range(10000):
        q = q_asimov(f, L, a, b, c)
        if q >= q_target or math.isinf(q):
            # bracket between f_prev (q < target) and f (q >= target or inf)
            lo, hi = (min(f_prev, f), max(f_prev, f))
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                qmid = q_asimov(mid, L, a, b, c)
                if qmid < q_target:
                    lo = mid
                else:
                    hi = mid
                if abs(hi - lo) < 1e-6:
                    break
            return 0.5 * (lo + hi)
        f_prev = f
        f *= 1.2  # gentle geometric growth to avoid overshooting
    raise RuntimeError(f"Did not find crossing on side {side}")

def asimov_limits_95(L, a, b, c):
    """
    Two-sided 95% CL limits (q=3.84) for the single-bin no-background case.
    Returns (f_minus, f_plus) in TeV^-4.
    """
    f_minus = _limit_one_side(L, a, b, c, side="-")
    f_plus  = _limit_one_side(L, a, b, c, side="+")
    return f_minus, f_plus

def fisher_approx_95(L, a, c):
    """
    Symmetric Fisher/Gaussian approximation near f=0 (rate-only, single bin):
      F = L * a^2 / c,  |f|_95 ≈ 1.96 / sqrt(F)
    Useful as a quick cross-check.
    """
    F = L * (a*a) / c
    return 1.96 / math.sqrt(F)

# ---- Example with your numbers ----
if __name__ == "__main__":
    # L in fb^-1; a,b,c in fb (f in TeV^-4)
    L = 1000.0
    a = -0.014647       # fb / (TeV^-4)
    b =  0.00089293     # fb / (TeV^-4)^2
    c = 15.3833         # fb

    f_lo, f_hi = asimov_limits_95(L, a, b, c)
    f95_fisher = fisher_approx_95(L, a, c)

    print(f"Two-sided 95% CL (Asimov): f in [{f_lo:.3f}, +{f_hi:.3f}] TeV^-4")
    print(f"Fisher symmetric approx:  |f|_95 ≈ {f95_fisher:.3f} TeV^-4")
