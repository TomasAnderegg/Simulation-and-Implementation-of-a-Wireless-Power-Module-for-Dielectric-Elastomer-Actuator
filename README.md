# ⚡ Wireless Power Module for Dielectric Elastomer Actuator

### Semester Project -- LAI Laboratory (EPFL)

**Author:** Tomas Garate Anderegg\
**Supervisor:** Maribel Caceres Rivera\
**Professors:** Yves Perriard & Paolo Germano

------------------------------------------------------------------------

# 📌 Project Overview

Design and experimental validation of a **13.56 MHz Wireless Power
Transfer (WPT) system** for powering a Dielectric Elastomer Actuator
(DEA) in biomedical applications.

**Target Power:** 300--500 mW\
**Max Distance:** 10 mm\
**Topology:** NRIC -- Serial-Serial Compensation

------------------------------------------------------------------------

# 🧠 Design Workflow

Simulation → Coil Design → PCB Layout → Fabrication → Measurement →
Efficiency Evaluation

------------------------------------------------------------------------

# 1️⃣ Resonant Circuit Parameters

  Parameter            Initial Value   Adjusted Value
  -------------------- --------------- ----------------
  Inductance (µH)      2.36            1.75
  Capacitance (pF)     53.3            73.3
  Quality Factor (Q)   201             149

Simulated in **LTSpice**.

------------------------------------------------------------------------

# 2️⃣ Final Coil Geometry

  Parameter         Value
  ----------------- ---------
  Inductance        1.84 µH
  Number of Turns   11
  Outer Diameter    40.4 mm
  Track Width       0.97 mm
  Spacing           0.87 mm

Designed using MATLAB + Current Sheet Approximation.

------------------------------------------------------------------------

# 3️⃣ Measured Electrical Results (Unloaded, 10 mm)

  Quantity                    Measured Value
  --------------------------- ----------------
  S11 (dB)                    -22.39
  Estimated Input Impedance   58.2 Ω
  Coupling Factor (K)         0.345

------------------------------------------------------------------------

# 4️⃣ Loaded Performance (1 kΩ Load, 10 mm)

  Quantity               Value
  ---------------------- -----------
  S21 (dB)               5.02
  S11 (dB)               -26.77
  Output Current (RMS)   8.91 mA
  Power Delivered        158.76 mW
  Efficiency             57.48 %

------------------------------------------------------------------------

# 📊 Theoretical vs Measured Comparison

  Parameter     Theoretical   Measured
  ------------- ------------- -----------
  Inductance    1.84 µH       1.67 µH
  Power (1kΩ)   300--500 mW   158.76 mW
  Efficiency    \~51 %        57.48 %

------------------------------------------------------------------------

# ⚠️ Key Limitations

-   Load not included in early resonance modeling\
-   Skin & proximity effects neglected\
-   Transmission line effects ignored\
-   Limited iterative tuning

------------------------------------------------------------------------

# 🛠 Tools Used

-   MATLAB\
-   LTSpice\
-   KiCAD\
-   Agilent 4294A Impedance Analyzer\
-   Agilent N5242A Network Analyzer

------------------------------------------------------------------------

# 👨‍🔬 Acknowledgment

The LAI Laboratory warmly welcomed me into this hardware-focused
project. I am extremely grateful for the opportunity, as it allowed me
to develop strong hands-on experience in RF design, PCB development, and
high-frequency measurements.

------------------------------------------------------------------------

# 📜 Conclusion

✔ Functional wireless power transfer at 13.56 MHz\
✔ Coupling factor K ≈ 0.345\
✔ Efficiency ≈ 57 %

This project highlighted the importance of load-aware design,
high-frequency modeling, and iterative tuning in WPT systems.
