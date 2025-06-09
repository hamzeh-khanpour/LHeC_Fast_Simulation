import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# =============================
# 1. Load the data
# =============================
filename = "cross_sections_FM2.txt"
fm2_vals = []
xsec_vals = []

with open(filename, 'r') as file:
    for line in file:
        if line.startswith("FM2") or not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                fm2 = float(parts[0].replace('e-12', '')) * 1e-12
                xsec = float(parts[1])
                fm2_vals.append(fm2)
                xsec_vals.append(xsec)
            except ValueError:
                continue

fm2_vals = np.array(fm2_vals)
xsec_vals = np.array(xsec_vals)
fm2_scaled = fm2_vals / 1e-12  # rescale for fitting



# =============================
# 2. Define fit function
# =============================
def quad_func(x, a, b):
    return a * x**2 + b



# =============================
# 3. Perform fit
# =============================
popt, pcov = curve_fit(quad_func, fm2_scaled, xsec_vals)
a, b = popt
perr = np.sqrt(np.diag(pcov))  # uncertainties



# =============================
# 4. Plot: Data + Fit + Residuals
# =============================
x_fit = np.linspace(-1100, 1100, 500)
y_fit = quad_func(x_fit, *popt)



# Plot main curve
plt.figure(figsize=(8, 5))
plt.plot(fm2_scaled, xsec_vals, 'o', label='Data')
plt.plot(x_fit, y_fit, '-', label='Fit: $a x^2 + b$')

fit_text = f"Fit coefficients:\n" \
           f"$a = ({a:.3e} \pm {perr[0]:.1e})$\n" \
           f"$b = ({b:.3e} \pm {perr[1]:.1e})$"
plt.text(0.05, 0.95, fit_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

plt.xlabel('FM2 Value (scaled by $10^{-12}$)')
plt.ylabel('Cross-section (pb)')
plt.title('Cross-section vs FM2 (Quadratic Fit)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("fit_cross_section_vs_FM2.pdf")
plt.show()



# Residuals
residuals = xsec_vals - quad_func(fm2_scaled, *popt)
plt.figure(figsize=(8, 3))
plt.axhline(0, color='gray', linestyle='--')
plt.plot(fm2_scaled, residuals, 'o')
plt.xlabel('FM2 Value (scaled by $10^{-12}$)')
plt.ylabel('Residuals (pb)')
plt.title('Residuals of the Fit')
plt.grid(True)
plt.tight_layout()
plt.savefig("fit_residuals.pdf")
plt.show()



# =============================
# 5. Save Fit Function
# =============================
def fitted_cross_section(fm2):
    """Compute cross-section from FM2 in e-12 units"""
    x = fm2 / 1e-12
    return quad_func(x, *popt)



# =============================
# 6. Prediction Table
# =============================
fm2_pred = np.array([-2000, -1000, -500, -100, 0, 100, 500, 1000, 2000]) * 1e-12
x_pred = fm2_pred / 1e-12
y_pred = quad_func(x_pred, *popt)



# Error propagation: dy = sqrt((∂f/∂a)^2 * σ_a^2 + (∂f/∂b)^2 * σ_b^2)
y_err = np.sqrt((x_pred**2)**2 * pcov[0, 0] + pcov[1, 1] + 2 * x_pred**2 * pcov[0, 1])


print("\nPrediction Table (including uncertainties):")
print(f"{'FM2 (e-12)':>12} {'σ (pb)':>12} {'± Δσ (pb)':>12}")
for x, y, err in zip(x_pred, y_pred, y_err):
    print(f"{x:12.1f} {y:12.5f} {err:12.5f}")



