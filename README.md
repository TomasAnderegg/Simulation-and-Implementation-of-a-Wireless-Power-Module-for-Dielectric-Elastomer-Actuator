# ⚡ Wireless Power Module for Dielectric Elastomer Actuator

### Semester Project — LAI Laboratory (EPFL)

**Section de Microtechnique — Spring 2024/2025**  
**Author:** Tomas Garate Anderegg  
**Supervisor:** Maribel Caceres Rivera  
**Professors:** Yves Perriard & Paolo Germano

---

## 📌 Project Overview

This semester project focused on the **simulation, design, and experimental validation of a Wireless Power Transfer (WPT) system** operating at **13.56 MHz** for biomedical applications.

The objective was to develop a **wireless implantable power module** capable of transferring energy from an external transmitter (Tx) to an implanted receiver (Rx) to power a **Dielectric Elastomer Actuator (DEA)**.

- **Target power delivery:** 300–500 mW  
- **Maximum separation distance:** ~10 mm  
- **Operating frequency:** 13.56 MHz

---

## 🧠 Technical Approach

### Selected Topology

After reviewing different WPT techniques:

- Non-Radiative Inductive Coupling (NRIC)
- Non-Radiative Magnetic Resonant Coupling (NRMRC)
- Acoustic Power Transfer (APT)

I selected:

> **NRIC with Serial-Serial (SS) compensation topology**

Reasons:
- Simpler architecture
- Well documented in literature
- Suitable for short-distance biomedical applications
- Good Power Transfer Efficiency (PTE)

---

## 🔄 Design Methodology

The project followed a structured hardware development pipeline:

### 1️⃣ Resonant Circuit Design

Target resonance frequency:

> f₀ = 13.56 MHz

Resonance condition:

> f₀ = 1 / (2π√LC)

Initial theoretical parameters:

| Parameter   | Initial  | Adjusted |
|-------------|----------|----------|
| Inductance  | 2.36 µH  | 1.75 µH  |
| Capacitance | 53.3 pF  | 73.3 pF  |
| Q Factor    | 201      | 149      |

Simulations performed in **LTSpice**.

---

### 2️⃣ Coil Geometry Design

A MATLAB script was developed to:

- Compute inductance using Current Sheet Approximation
- Sweep geometric parameters:
  - Number of turns (N)
  - Track width (w)
  - Inner/outer diameter

Final coil parameters:

| Parameter      | Value    |
|----------------|----------|
| Inductance     | 1.84 µH  |
| Turns (N)      | 11       |
| Outer Diameter | 40.4 mm  |
| Track Width    | 0.97 mm  |
| Spacing        | 0.87 mm  |

---

### 3️⃣ PCB Design

- Designed using **KiCAD**
- Included test pads for tuning
- Connector outputs for measurement
- Matching topology flexibility

---

### 4️⃣ Experimental Validation

#### Inductance Measurement

Measured with **Agilent 4294A Precision Impedance Analyzer**.

Measured inductance at 13.56 MHz:

> L ≈ 1.67 µH — Relative error ≈ 9%

Measured coil resistance:

> R ≈ 1 Ω

---

### 5️⃣ Network Analysis

Measured using **Agilent N5242A PNA-X Network Analyzer**.

Measured S-parameters (10 mm separation, unloaded):

> S11 = -22.39 dB  
> Estimated input impedance ≈ 58.2 Ω

Estimated coupling factor:

> K ≈ 0.345

---

## 🔌 Loaded System Performance (1 kΩ Load)

At 13.56 MHz and 10 mm separation:

| Parameter   | Measured       |
|-------------|----------------|
| I₂ (RMS)    | 8.91 mA        |
| PLoad       | 158.76 mW      |
| Efficiency  | **57.48%**     |

---

## 📊 Comparative Results

| Parameter    | Theoretical  | Measured       |
|--------------|--------------|----------------|
| Inductance   | 1.84 µH      | 1.67 µH        |
| PLoad (1kΩ)  | 300–500 mW   | 158.76 mW      |
| Efficiency   | ~51%         | **57.48%**     |

Some parameters exceeded expectations, while delivered power fell below the initial target range.

---

## ⚠️ Key Limitations Identified

- Skin depth & proximity effects neglected
- Load not included in early resonance design
- Transmission line effects ignored
- Limited fine-tuning iterations
- High-frequency parasitic effects underestimated

These factors contributed to deviations between theoretical and measured performance.

---

## 🧪 Efficiency Modeling

Using measured S-parameters, efficiency was computed using two-port network theory:

> η = PLoad / Pin

For K ≈ 0.3, measured efficiency:

> **η ≈ 57.48%**

Which is consistent with LTSpice simulations of maximum achievable efficiency vs coupling factor.

---

## 💡 Lessons Learned

- Load integration must be included early in design
- High-frequency parasitics significantly impact performance
- WPT system tuning is inherently iterative
- PCB geometry and component tolerances matter greatly at MHz frequencies

---

## 🚀 Future Improvements

- Include load in resonance modeling from beginning
- Account for skin and proximity effects
- Improve impedance matching
- Perform iterative tuning cycles
- Refine coil optimization for better PTE

---

## 🛠 Tools & Technologies

- MATLAB
- LTSpice
- KiCAD
- Agilent 4294A Impedance Analyzer
- Agilent N5242A Network Analyzer

---

## 👨‍🔬 Acknowledgment

I wanted to work on a hardware-focused project, and the **LAI Laboratory** welcomed me warmly. I am truly grateful for this opportunity — it allowed me to gain substantial hands-on experience in:

- RF design
- PCB development
- Resonant circuits
- S-parameter analysis
- Biomedical WPT systems

This project significantly strengthened both my theoretical understanding and practical engineering skills.

---

## 📜 Conclusion

This project demonstrated the full pipeline of a real-world hardware development process:

> Simulation → Coil Design → PCB Fabrication → Measurement → Network Analysis → Efficiency Evaluation

While theoretical assumptions provided a strong starting point, experimental validation revealed the importance of load-aware design, high-frequency modeling, and iterative tuning.

**Final system achieved:**

✔ Measured coupling factor K ≈ 0.345  
✔ Efficiency ≈ 57%  
✔ Functional wireless power transfer at 13.56 MHz

This work lays the foundation for further optimization toward high-efficiency biomedical WPT systems.
