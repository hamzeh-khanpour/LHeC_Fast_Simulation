import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data from file
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

# Convert to arrays
fm2_vals = np.array(fm2_vals)
xsec_vals = np.array(xsec_vals) * 1000.0
fm2_scaled = fm2_vals / 1e-12  # rescale

# Fit function
def symmetric_poly(x, a, b, c):
    return a * x + b * x**2 + c

# Perform fit
popt, pcov = curve_fit(symmetric_poly, fm2_scaled, xsec_vals)
a, b, c = popt

# Create fit curve
x_fit = np.linspace(-1100, 1100, 500)
y_fit = symmetric_poly(x_fit, *popt)

# Plot
plt.figure(figsize=(8, 6))
plt.plot(fm2_scaled, xsec_vals, 'o', label='Data')
plt.plot(x_fit, y_fit, '-', label='Fit: $a x + b x^2 + c$')

# Add fitted values as text box
fit_text = f"Fit coefficients:\n" \
           f"$a = {a:.3e}$\n" \
           f"$b = {b:.3e}$\n" \
           f"$c = {c:.3e}$"
plt.text(0.05, 0.95, fit_text, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

# Labels and style
plt.xlabel('$f_{M2}/\Lambda^4$ Value (scaled by $10^{-12}$)')
plt.ylabel('Cross-section (fb)')
plt.title('Cross-section vs $f_{M2}/\Lambda^4$ (Fit)')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save plot as PDF
plt.savefig("fit_cross_section_vs_FM2.pdf", format='pdf')

# Show plot
plt.show()

# Print fitted parameters to console
print("Fitted parameters:")
print(f"a = {a:.6e}")
print(f"b = {b:.6e}")
print(f"c = {c:.6e}")
