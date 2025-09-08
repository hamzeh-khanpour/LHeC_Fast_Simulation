import math

# -----------------------------
# Inputs (fiducial, in pb)
# -----------------------------
data = {
    400: {"c": 1.0075e-02, "sp": 1.1160e-02, "sm": 1.1377e-02},
    500: {"c": 3.5750e-03, "sp": 4.3741e-03, "sm": 4.4594e-03},
    600: {"c": 1.2235e-03, "sp": 1.6982e-03, "sm": 1.7313e-03},
    700: {"c": 4.1919e-04, "sp": 6.3510e-04, "sm": 6.4748e-04},
    800: {"c": 1.3810e-04, "sp": 2.1418e-04, "sm": 2.1836e-04},
}

# Analysis settings
L_fb = 1000.0          # integrated luminosity in fb^-1
EPS = 1.0              # efficiency (set to 0.20 to include your 20% later if desired)
PB_TO_FB = 1e3         # conversion: 1 pb = 1000 fb
CL95_Q = 3.84          # 95% CL for 1 d.o.f.

# -----------------------------
# EFT coefficient extraction
# -----------------------------
def extract_AB(c, sp, sm, f0=10.0):
    """
    Given c=σ(0), sp=σ(+f0), sm=σ(-f0) in pb and f0 in TeV^-4,
    return (A, B) in pb/TeV^-4 and pb/TeV^-8 respectively.
    """
    A = (sp - sm) / (2.0 * f0)
    B = (sp + sm - 2.0 * c) / (2.0 * f0 * f0)
    return A, B

# -----------------------------
# Gaussian/Asimov analytic tolerance
# -----------------------------
def sigma_tolerance_gaussian(c_pb, L_fb, eps):
    """
    S such that |σ(f) - c| <= S approximates the 95% CL (Gaussian).
    Units: pb.
    S = 1.96 * sqrt( c / (L * eps * 1000) )
    """
    K = L_fb * eps * PB_TO_FB
    return 1.96 * math.sqrt(c_pb / K)

def interval_gaussian(c, A, B, L_fb, eps):
    """
    Solve B f^2 + A f - S <= 0 for f, with S from Gaussian tolerance.
    Returns (f_minus, f_plus). Handles linear case B≈0.
    """
    S = sigma_tolerance_gaussian(c, L_fb, eps)
    if abs(B) < 1e-16:
        # purely linear
        if abs(A) < 1e-16:
            return (-float("inf"), float("inf"))  # flat case (should not happen here)
        f0 = S / abs(A)
        return (-f0, +f0) if A > 0 else (-f0, +f0)
    disc = A*A + 4.0*B*S
    if disc < 0:
        return (float("nan"), float("nan"))
    sqrt_disc = math.sqrt(disc)
    f_minus = (-A - sqrt_disc) / (2.0 * B)
    f_plus  = (-A + sqrt_disc) / (2.0 * B)
    # Ensure ordering
    if f_minus > f_plus:
        f_minus, f_plus = f_plus, f_minus
    return (f_minus, f_plus)

# -----------------------------
# Exact Poisson Asimov solution
# -----------------------------
def mu_from_sigma(sigma_pb, L_fb, eps):
    """Expected counts μ = L * eps * σ * 1000 (since 1 pb = 1000 fb)."""
    return L_fb * eps * sigma_pb * PB_TO_FB

def q_asimov(mu, mu0):
    """One-bin Asimov test statistic."""
    if mu <= 0 or mu0 <= 0:
        return float("inf")
    return 2.0 * (mu - mu0 + mu0 * math.log(mu0 / mu))

def interval_asimov_exact(c, A, B, L_fb, eps, q_target=CL95_Q):
    """
    Solve q(f) = q_target for the two roots around f=0 using bracketing+bisection.
    Returns (f_minus, f_plus).
    """
    def qf(f):
        sigma = c + A*f + B*f*f
        return q_asimov(mu_from_sigma(sigma, L_fb, eps), mu_from_sigma(c, L_fb, eps))

    # Find right root
    f_lo, f_hi = 0.0, 0.0
    q_lo = qf(f_lo)
    step = 0.1
    # march outwards to the right until we bracket
    while True:
        f_hi = f_lo + step
        q_hi = qf(f_hi)
        if q_hi >= q_target:
            break
        step *= 1.8
        if step > 1e6:  # safety
            return (float("nan"), float("nan"))
    # bisection on the right
    for _ in range(80):
        f_mid = 0.5*(f_lo + f_hi)
        q_mid = qf(f_mid)
        if q_mid < q_target:
            f_lo = f_mid
        else:
            f_hi = f_mid
    f_plus = 0.5*(f_lo + f_hi)

    # Find left root (mirror search)
    f_hi, f_lo = 0.0, 0.0
    q_hi = qf(f_hi)
    step = -0.1
    while True:
        f_lo = f_hi + step
        q_lo = qf(f_lo)
        if q_lo >= q_target:
            break
        step *= 1.8
        if abs(step) > 1e6:
            return (float("nan"), float("nan"))
    # bisection on the left
    for _ in range(80):
        f_mid = 0.5*(f_lo + f_hi)
        q_mid = qf(f_mid)
        if q_mid < q_target:
            f_hi = f_mid
        else:
            f_lo = f_mid
    f_minus = 0.5*(f_lo + f_hi)

    return (f_minus, f_plus)

# -----------------------------
# Run and print tables
# -----------------------------
print(f"Using L = {L_fb:.0f} fb^-1, efficiency ε = {EPS:.2f} (set EPS=1.0 for 'without 20%').\n")

print("EFT coefficients per W_min (c in pb, A in pb/TeV^-4, B in pb/TeV^-8):")
print(f"{'W_min[GeV]':>10}  {'c':>12}  {'A':>12}  {'B':>12}")
AB = {}
for W, d in data.items():
    c, sp, sm = d["c"], d["sp"], d["sm"]
    A, B = extract_AB(c, sp, sm, f0=10.0)
    AB[W] = (c, A, B)
    print(f"{W:>10}  {c:12.6e}  {A:12.6e}  {B:12.6e}")

print("\n95% CL intervals for f (TeV^-4) — Gaussian (analytic) and exact Asimov:")
print(f"{'W_min[GeV]':>10}  {'Gaussian [f-, f+]':>30}  {'Asimov-exact [f-, f+]':>32}")
for W in sorted(AB):
    c, A, B = AB[W]
    fG = interval_gaussian(c, A, B, L_fb, EPS)
    fE = interval_asimov_exact(c, A, B, L_fb, EPS, q_target=CL95_Q)
    print(f"{W:>10}  [{fG[0]:7.2f}, {fG[1]:7.2f}]        [{fE[0]:7.2f}, {fE[1]:7.2f}]")

# Optional: show the interval center and width
print("\nCenters f* = -A/(2B) (TeV^-4) for intuition:")
for W, (c, A, B) in AB.items():
    fstar = -A/(2.0*B) if abs(B) > 0 else float('nan')
    print(f"W>{W:>3} GeV : f* ≈ {fstar:.3f}")
