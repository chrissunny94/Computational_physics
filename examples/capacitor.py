#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from rk import RK4

# ------------------------------------
# Circuit parameters
# ------------------------------------
R = 1e3          # Resistance (ohm)
C = 40e-3        # Capacitance (F)
C0 = 0           # Initial charge
U0 = 5e3         # Voltage source (V)

# dQ/dt = U/R – Q/(RC)
Qdot = lambda t, Q: U0/R - Q/(R*C)

# ------------------------------------
# Solve using RK4
# ------------------------------------
rk = RK4(Qdot)
t, y = rk.solve([C0], 0.1, 500)
Q = y[0]

# ------------------------------------
# STATIC VISUALIZATION
# ------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(t, Q, label="Capacitor Charge Q(t)")
plt.title("Capacitor Charging Visualization")
plt.xlabel("Time (s)")
plt.ylabel("Charge (C)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
