# to calculate limits for 1000 fb⁻¹ and assuming ‘observed’ x-secs and no backgrounds
# Luminosity scenario: Projected sensitivity at an integrated luminosity of 𝐿=1000 fb−1 (FCChe-scale).
# Observed cross section (“Asimov data”): Take the SM prediction as the observed measurement
# No backgrounds: Pretend your measurement is background-free. That makes the statistics purely Poisson on the signal yield and gives an optimistic / best-case limit.

import math

def q_asimov(f, L, a, b, c):
    mu0 = L * c
    muf = L * (c + a*f + b*f*f)
    if muf <= 0:
        return float('inf')
    return 2.0 * (muf - mu0 + mu0 * math.log(mu0 / muf))

def _limit_one_side(L, a, b, c, q_target=3.84, side="+"):
    # start from f=0 and march outward
    step = 1.0 if side == "+" else -1.0
    f_prev = 0.0
    f = step
    for _ in range(10000):
        q = q_asimov(f, L, a, b, c)
        if q >= q_target or math.isinf(q):
            lo, hi = (min(f_prev, f), max(f_prev, f))
            # --- minimal fix: preserve bracket & tighten tolerance ---
            qlo = q_asimov(lo, L, a, b, c)
            for _ in range(500):  # was 200
                mid = 0.5 * (lo + hi)
                qmid = q_asimov(mid, L, a, b, c)
                if (qmid - q_target) * (qlo - q_target) > 0:
                    lo, qlo = mid, qmid
                else:
                    hi = mid
                if abs(hi - lo) < 1e-10:  # was 1e-6
                    break
            return 0.5 * (lo + hi)
        f_prev = f
        f *= 1.2
    raise RuntimeError(f"Did not find crossing on side {side}")

def asimov_limits_95(L, a, b, c):
    f_minus = _limit_one_side(L, a, b, c, side="-")
    f_plus  = _limit_one_side(L, a, b, c, side="+")
    return f_minus, f_plus

def fisher_approx_95(L, a, c):
    F = L * (a*a) / c
    return 1.96 / math.sqrt(F)

# ---- Example with your numbers ----
if __name__ == "__main__":
    L =  1000.0
    a = -1.469751e-01         # fb / (TeV^-4)
    b =  7.055479e-01         # fb / (TeV^-4)^2
    c =  1.083117e+02         # fb

    f_lo, f_hi = asimov_limits_95(L, a, b, c)
    f95_fisher = fisher_approx_95(L, a, c)

    print(f"Two-sided 95% CL (Asimov): f in [{f_lo:.6f}, +{f_hi:.6f}] TeV^-4")
    print(f"Fisher symmetric approx:  |f|_95 ≈ {f95_fisher:.3f} TeV^-4")
