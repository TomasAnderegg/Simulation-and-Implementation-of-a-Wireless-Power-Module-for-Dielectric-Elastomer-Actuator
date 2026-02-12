# -*- coding: utf-8 -*-
"""
S-Parameter Analysis at 13.56 MHz: Current, Power, and Z21 Computation
@author: adapted
"""

import skrf as rf
import matplotlib.pyplot as plt
import numpy as np
from skrf.plotting import plot_smith

# === Load .s2p file ===
network = rf.Network("tr.s2p")  # Replace with your actual file

# === Parameters ===
Z0 = 50
ZL = 50 + 0j
a1 = 1
f = network.f
f_MHz = f / 1e6

# === Extract S-parameters ===
S11 = network.s[:, 0, 0]
S21 = network.s[:, 1, 0]
S12 = network.s[:, 0, 1]
S22 = network.s[:, 1, 1]

# === Load reflection coefficient ===
Gamma_L = (ZL - Z0) / (ZL + Z0)

# === Reflected wave at port 2 ===
b2 = (S21 * a1) / (1 - S22 * Gamma_L)
a2 = Gamma_L * b2

# === Current into the load ===
I2 = (a2 - b2) / np.sqrt(Z0)
I2_rms = np.abs(I2) / np.sqrt(2)

# === Power delivered to the load ===
P_L = np.abs(b2)**2 * (1 - np.abs(Gamma_L)**2)

# === Loaded S21 definition ===
S21_loaded = S21 * (1 + Gamma_L) / (1 - S22 * Gamma_L)

# === Input reflection coefficient and power ===
Gamma_in = S11 + (S12 * S21 * Gamma_L) / (1 - S22 * Gamma_L)
Z_in = Z0 * (1 + Gamma_in) / (1 - Gamma_in)
P_in = (np.abs(a1)**2) * (1 - np.abs(Gamma_in)**2)

# === Efficiency ===
eta = (P_L / P_in) * 100

# === Transform S -> Z matrix ===
Z_matrices = []
I = np.identity(2)

for i in range(len(f)):
    S = np.array([[S11[i], S12[i]], [S21[i], S22[i]]])
    Z = Z0 * (I + S) @ np.linalg.inv(I - S)
    Z_matrices.append(Z)

Z_matrices = np.array(Z_matrices)
Z21 = Z_matrices[:, 1, 0]

# === Results at 13.56 MHz ===
f_target = 13.56e6
idx = np.argmin(np.abs(f - f_target))

print(f"--- Results at {f[idx]/1e6:.3f} MHz ---")
print(f"Z_in (input impedance) = {np.abs(Z_in[idx]):.2f} Ω ∠ {np.angle(Z_in[idx], deg=True):.2f}°")
print(f"|S21_loaded|         = {20*np.log10(np.abs(S21_loaded[idx])):.2f} dB")
print(f"|S11|                = {20*np.log10(np.abs(S11[idx])):.2f} dB")
print(f"Γ_L (magnitude)      = {np.abs(Gamma_L):.4f} ∠ {np.angle(Gamma_L, deg=True):.2f}°")
print(f"I2 (RMS)             = {I2_rms[idx]*1000:.2f} mA")
print(f"P_L (delivered)      = {P_L[idx]*1000:.2f} mW")
print(f"Efficiency η         = {eta[idx]:.2f} %")
print(f"Z21 = {np.abs(Z21[idx]):.2f} Ω ∠ {np.angle(Z21[idx], deg=True):.2f}°")
M12 = (np.abs(Z21[idx]) * np.abs(np.sin(np.angle(Z21[idx]))) / (2 * np.pi * f_target)) # Mutual inductance
print(f"M12 = {M12:.2e}")
K = M12/1.67053e-6 # Coupling Factor
print(f"K = {K:.2e}")

# === Plots ===
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(f_MHz, 20*np.log10(np.abs(S21_loaded)), label='|S21_loaded| (dB)')
plt.plot(f_MHz, 20*np.log10(np.abs(S11)), label='|S11| (dB)', linestyle='--')
plt.ylabel('Magnitude (dB)')
plt.title(f'Loaded Response with ZL = {ZL} Ω')
plt.legend(); plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(f_MHz, I2_rms*1000, color='red')
plt.ylabel('Current I2 (mA RMS)')
plt.title('RMS Current into Load')
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(f_MHz, P_L*1000, color='green')
plt.xlabel('Frequency (MHz)')
plt.ylabel('Delivered Power (mW)')
plt.title('Power Delivered to Load')
plt.grid(True)

plt.tight_layout()
plt.show()
