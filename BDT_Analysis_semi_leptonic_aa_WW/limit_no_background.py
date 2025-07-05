# ================================================================================
#   Hamzeh Khanpour — June 2025
#   No-Background 95% CL Limit on FM2 (cut-and-count + profile likelihood)
# ================================================================================

import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# --- Inputs ---
luminosity_fb = 100.0
delta_sys = 0.10  # 10% systematic uncertainty

signal_efficiency = 0.1405  # Example from your latest BDT output
background_efficiency  = 0.1205

sigma_sm_fb = 1.538327e+01 * background_efficiency # SM cross section in fb

# --- Cross-section function (in fb) ---
def sigma_fm2_fb(fm2):
    a = -1.464725e-02
    b = 8.929319e-04
    c = 1.538327e+01
    return a * fm2 + b * fm2**2 + c

# --- Derived ---
n_sm = luminosity_fb * sigma_sm_fb
delta_stat = 1.0 / np.sqrt(n_sm)
delta_tot = np.sqrt(delta_stat**2 + delta_sys**2)

# --- χ² test function ---
def chi2_stat(fm2):
    sigma_bsm = sigma_fm2_fb(fm2) * signal_efficiency
    return ((sigma_sm_fb - sigma_bsm) / (sigma_sm_fb * delta_tot)) ** 2

# --- Find limits where χ² = 3.84 ---
def chi2_diff(fm2):
    return chi2_stat(fm2) - 3.84

# Scan around fm2 = 0
result_lower = root_scalar(chi2_diff, bracket=[-100, 0.0001], method='brentq')
result_upper = root_scalar(chi2_diff, bracket=[0.0001, 100], method='brentq')

# --- Output ---
print("========= χ² Limit Setting at 95% CL =========")
print(f"Luminosity: {luminosity_fb} fb⁻¹")
print(f"SM cross section: {sigma_sm_fb:.4f} fb")
print(f"Stat. error: {delta_stat:.4f}, Total rel. error: {delta_tot:.4f}")
print(f"95% CL Exclusion Range for fM2: {result_lower.root:.3f} to {result_upper.root:.3f} [TeV⁻⁴]")

# --- Optional: Plot χ² curve ---
fm2_vals = np.linspace(-100, 100, 200)
chi2_vals = [chi2_stat(fm2) for fm2 in fm2_vals]

plt.plot(fm2_vals, chi2_vals, label=r"$\chi^2(f_{M2})$")
plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold')
plt.xlabel(r"$f_{M2} / \Lambda^4\ [\mathrm{TeV}^{-4}]$")
plt.ylabel(r"$\chi^2$")
plt.legend()
plt.grid(True)
plt.title(r"Profile $\chi^2$ for $f_{M2}$ at 95% CL")
plt.show()


# ================================================================================



import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# ====== Inputs ======
luminosity_fb = 100.0               # fb^-1
sigma_sm_fb = 15.38327               # SM cross section [fb]
delta_sys = 0.10                    # 10% systematic uncertainty
confidence_level = 0.95             # 95% CL
chi2_threshold = 3.84               # For 1 d.o.f. at 95% CL

signal_efficiency = 0.1405  # Example from your latest BDT output
background_efficiency  = 0.1205

# ====== Cross-section function from fit ======
def sigma_fm2_fb(fm2):
    a = -1.464725e-02
    b = 8.929319e-04
    c = 15.38327
    return a * fm2 + b * fm2**2 + c

# ====== Derived quantities ======
n_sm = luminosity_fb * sigma_sm_fb
delta_stat = 1 / np.sqrt(n_sm)
delta_tot = np.sqrt(delta_stat**2 + delta_sys**2)


# ====== χ² function ======
def chi2_fm2(fm2):
    sigma_bsm = sigma_fm2_fb(fm2)
    return ((sigma_sm_fb - sigma_bsm) / (sigma_sm_fb * delta_tot))**2


# ====== Root-finding to find limits ======
def chi2_diff(fm2):
    return chi2_fm2(fm2) - chi2_threshold


lower = root_scalar(chi2_diff, bracket=[-100, 0], method='brentq').root
upper = root_scalar(chi2_diff, bracket=[0, 100], method='brentq').root


# ====== Print result ======
print("=========== χ²-Based Limit (Paper Method) ===========")
print(f"√s = 1.98 TeV, Luminosity = {luminosity_fb} fb⁻¹")
print(f"σ_SM = {sigma_sm_fb:.4f} fb")
print(f"Total Relative Error: δ_tot = {delta_tot:.5f}")
print(f"95% CL Exclusion Range for fM2: {lower:.3f} to {upper:.3f} [TeV⁻⁴]")


# ====== Optional Plot ======
fm2_vals = np.linspace(-100, 100, 500)
chi2_vals = [chi2_fm2(fm2) for fm2 in fm2_vals]


plt.plot(fm2_vals, chi2_vals, label=r"$\chi^2(f_{M2})$")
plt.axhline(chi2_threshold, color='red', linestyle='--', label='95% CL threshold')
plt.xlabel(r"$f_{M2}/\Lambda^4\ [\mathrm{TeV}^{-4}]$")
plt.ylabel(r"$\chi^2$")
plt.title(r"Exclusion Scan using $\chi^2$ Method")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()




# ================================================================================
#   LHeC Semi-Leptonic WW → Limit on fM2 with 1000 fb⁻¹ (No-Background Scenario)
#   Includes:
#     - Chi² Method (as used in arXiv:2002.09748)
#     - Bayesian Method (Poisson 95% CL with 0 background)
#   Author: Hamzeh Khanpour — Updated July 2025
# ================================================================================

import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# ========================
# 🔢 Input Parameters
# ========================

luminosity_fb = 100.0               # Integrated luminosity [fb⁻¹]
eff_signal = 0.15                    # Total signal efficiency (acceptance × selection × BDT)
n_limit_bayes = 3.0                  # 95% CL Bayesian limit for 0 observed events
delta_sys = 0.10                     # 10% systematic uncertainty
sigma_sm_fb = 15.38327               # SM cross-section [fb] for aa → WW (after acceptance cuts)

# Cross-section function σ(fM2) [fb]
def sigma_fm2_fb(fm2):
    a = -1.464725e-02
    b = 8.929319e-04
    c = 15.38327
    return a * fm2 + b * fm2**2 + c


# ========================
# 📊 Chi² Method
# ========================

# Observed signal yield under SM
n_sm = sigma_sm_fb * luminosity_fb

# Combined uncertainty
delta_stat = 1.0 / np.sqrt(n_sm)                 # statistical
delta_tot = np.sqrt(delta_stat**2 + delta_sys**2)

def chi2_stat(fm2):
    sigma_bsm = sigma_fm2_fb(fm2)
    return ((sigma_sm_fb - sigma_bsm) / (sigma_sm_fb * delta_tot))**2

def chi2_diff(fm2):
    return chi2_stat(fm2) - 3.84

try:
    fm2_lower_chi2 = root_scalar(chi2_diff, bracket=[-100, -1e-3], method='brentq').root
    fm2_upper_chi2 = root_scalar(chi2_diff, bracket=[1e-3, 100], method='brentq').root
    chi2_success = True
except:
    chi2_success = False


# ========================
# 📏 Bayesian No-Background Limit
# ========================

# Bayesian cross-section limit
sigma_limit_fb = n_limit_bayes / (luminosity_fb * eff_signal)

def diff_bayes(fm2):
    return sigma_fm2_fb(fm2) - sigma_limit_fb

try:
    fm2_lower_bayes = root_scalar(diff_bayes, bracket=[-100, -1e-3], method='brentq').root
    fm2_upper_bayes = root_scalar(diff_bayes, bracket=[1e-3, 100], method='brentq').root
    bayes_success = True
except:
    bayes_success = False


# ========================
# 📈 Plot χ² Curve (Optional)
# ========================
plot_chi2 = True
if plot_chi2:
    fm2_vals = np.linspace(-100, 100, 500)
    chi2_vals = [chi2_stat(fm2) for fm2 in fm2_vals]

    plt.plot(fm2_vals, chi2_vals, label=r"$\chi^2(f_{M2})$")
    plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold')
    if chi2_success:
        plt.axvline(fm2_lower_chi2, color='gray', linestyle=':', label=f"Lower: {fm2_lower_chi2:.2f}")
        plt.axvline(fm2_upper_chi2, color='gray', linestyle=':', label=f"Upper: {fm2_upper_chi2:.2f}")
    plt.xlabel(r"$f_{M2}/\Lambda^4$ [$\mathrm{TeV}^{-4}$]")
    plt.ylabel(r"$\chi^2$")
    plt.grid(True)
    plt.legend()
    plt.title("Profile χ² Scan for $f_{M2}$ at 95% CL")
    plt.tight_layout()
    plt.savefig("limit_chi2_scan_fm2.pdf")
    plt.show()



# ========================
# 🖨️ Final Report
# ========================
print("=========== LIMIT RESULTS (1000 fb⁻¹) ===========")
print(f"Signal efficiency: {eff_signal:.3f}")
print(f"SM Cross section (σ_SM): {sigma_sm_fb:.4f} fb")
print(f"Stat error: {delta_stat:.4f}, Total relative error: {delta_tot:.4f}")
print(f"95% CL Bayesian limit (σ_limit): {sigma_limit_fb:.4f} fb")

# χ² method result
if chi2_success:
    print(f"🟦 χ² Excluded FM2 Range: [{fm2_lower_chi2:.3f}, {fm2_upper_chi2:.3f}] TeV⁻⁴")
else:
    print("❌ χ² root-finding failed — no intersection found.")

# Bayesian method result
if bayes_success:
    print(f"🟨 Bayesian No-Background FM2 Range: [{fm2_lower_bayes:.3f}, {fm2_upper_bayes:.3f}] TeV⁻⁴")
else:
    print("❌ Bayesian root-finding failed — σ(fM2) does not reach the 95% CL limit.")




