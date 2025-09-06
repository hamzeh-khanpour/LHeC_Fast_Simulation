import math
import numpy as np
import matplotlib.pyplot as plt

def q_asimov(f, L, a, b, c):
    mu0 = L * c
    muf = L * (c + a*f + b*f*f)
    if muf <= 0:
        return float('inf')
    return 2.0 * (muf - mu0 + mu0 * math.log(mu0 / muf))

def _limit_one_side(L, a, b, c, q_target=3.84, side="+"):
    step = 1.0 if side == "+" else -1.0
    f_prev = 0.0
    f = step
    for _ in range(10000):
        q = q_asimov(f, L, a, b, c)
        if q >= q_target or math.isinf(q):
            lo, hi = (min(f_prev, f), max(f_prev, f))
            for _ in range(200):
                mid = 0.5*(lo+hi)
                qmid = q_asimov(mid, L, a, b, c)
                if qmid < q_target:
                    lo = mid
                else:
                    hi = mid
                if abs(hi-lo) < 1e-6:
                    break
            return 0.5*(lo+hi)
        f_prev = f
        f *= 1.2
    raise RuntimeError("Did not find crossing")

def asimov_limits_95(L, a, b, c):
    return _limit_one_side(L,a,b,c,side="-"), _limit_one_side(L,a,b,c,side="+")

# --- Physics parameters (your fit results) ---
a = -0.014647     # fb / (TeV^-4)
b =  0.00089293   # fb / (TeV^-4)^2
c = 15.3833       # fb

# --- Scan over luminosity ---
L_values = np.logspace(1, 4, 40)  # from 10 to 10000 fb^-1
limits_minus = []
limits_plus  = []

for L in L_values:
    f_lo, f_hi = asimov_limits_95(L, a, b, c)
    limits_minus.append(f_lo)
    limits_plus.append(f_hi)

# --- Plot ---
plt.figure(figsize=(7,5))
plt.plot(L_values, limits_minus, label=r"$f_{M2}^-$ (95% CL)", color="blue")
plt.plot(L_values, limits_plus,  label=r"$f_{M2}^+$ (95% CL)", color="red")
plt.xscale("log")
plt.yscale("linear")
plt.xlabel("Integrated luminosity [fb$^{-1}$]", fontsize=12)
plt.ylabel(r"95% CL on $f_{M2}$ [TeV$^{-4}$]", fontsize=12)
plt.title("Expected 95% CL limits vs luminosity (Asimov, 1 bin, no backgrounds)")
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.legend()
plt.tight_layout()

# Save plot as PDF
plt.savefig("no-no-bkgs-plots.pdf", format='pdf')

plt.show()
