import numpy as np

def calculate_ceB(delta_a_tau, m_tau, v, e, Lambda, cos_theta_W):
    """
    Calculate the Wilson coefficient c_tau_B (ceB) for the anomalous magnetic moment of tau.

    Parameters:
    delta_a_tau (float): Anomalous magnetic moment of tau
    m_tau (float): Mass of tau lepton in GeV
    v (float): Vacuum expectation value of Higgs field in GeV
    e (float): Electromagnetic coupling constant in GeV
    Lambda (float): New physics scale in GeV
    cos_theta_W (float): Cosine of the weak mixing angle

    Returns:
    float: The value of c_tau_B (ceB)
    """
    return (delta_a_tau * e * Lambda**2) / (2 * m_tau * np.sqrt(2) * v * cos_theta_W)

# Given parameters
delta_a_tau = 0.0042  # Anomalous magnetic moment of tau
m_tau = 1.77686       # Mass of tau in GeV
v = 246               # Higgs VEV in GeV
e = 0.313             # Electromagnetic coupling constant in GeV
Lambda = 2000         # New physics scale in GeV (2 TeV)
cos_theta_W = 0.876   # Cosine of the weak mixing angle

# Compute ceB
ceB = calculate_ceB(delta_a_tau, m_tau, v, e, Lambda, cos_theta_W)
print(f"Calculated ceB: {ceB:.4f}")
