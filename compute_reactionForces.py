import pickle
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

with open("savedResults.pkl", 'rb') as handle:
    saved_results = pickle.load(handle)

t = saved_results["t_eval"]
y = saved_results["integratedCoordinates"]
LList = saved_results["lambdas"]
R_LG1List = saved_results["reactionForces"]
XddList = saved_results["Xdd"]

# Solved lambdas
L = np.zeros((4, len(LList)))
for i in range(len(LList)):
    L[:, i] = LList[i]

# Solved accelerations
Xdd = np.zeros((6, len(XddList)))
for i in range(len(XddList)):
    Xdd[:, i] = XddList[i]

# Computed reaction forces from Lagrange 1
R_LG1 = np.zeros((6, len(R_LG1List)))
for i in range(len(R_LG1List)):
    R_LG1[:, i] = R_LG1List[i]

# Geometry and inertia
m1 = 1           # kg
m2 = 1           # kg
g = 9.81         # m/s^2
l1 = 0.3         # m
l2 = 0.3         # m
J1 = m1*l1*l1/12 # kgm^2
J2 = m2*l2*l2/12 # kgm^2

# CoM accelerations
M = m1+m2
# x-Coordinate
xSdd = (m1*Xdd[0, :] + m2*Xdd[3, :]) / M
# y-Coordinate
ySdd = (m1*Xdd[1, :] + m2*Xdd[4, :]) / M

# Reaction forces from CoM dynamics
Fx = M*xSdd
Fy = M*ySdd + M*g

# Reaction forces from Lagrangian Equations Type I
Fx_LG1 = -(R_LG1[0, :] + R_LG1[3, :])
Fy_LG1 = -(R_LG1[1, :] + R_LG1[4, :])

# Plotting reactions forces from LG1 and CoM dynamics
title_text = f"Parameters: Pendulum 1 ({m1} kg, {l1} m) | Pendulum 2 ({m2} kg, {l2} m)"

tspan = np.linspace(t[0], t[-1], len(XddList))

fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.8), sharey=True)
axs[0].plot(tspan, Fx, label=r"$F_x$ (CoM)")
axs[0].plot(tspan, Fy, label=r"$F_y$ (CoM)")
axs[0].set_title(r"Reaction Forces with CoM Dynamics", loc="left")

axs[1].plot(tspan, Fx_LG1, label=r"$F_x$ (LG1)")
axs[1].plot(tspan, Fy_LG1, label=r"$F_y$ (LG1)")
axs[1].text(
    0.0, 105, r"$R=-{J_F}^T \lambda$", fontsize='large',
    backgroundcolor='white'
)
axs[1].set_title(r"Reaction Forces with Lagrangian Type I", loc="left")

axs[0].set_ylabel(r"$N$")
axs[0].grid(visible=True, which="major")
axs[1].grid(visible=True, which="major")

axs[0].set_xlabel(r"Time ($s$)")
axs[1].set_xlabel(r"Time ($s$)")
fig.suptitle(title_text, x=0.085, y=0.95, ha="left")

axs[0].legend()
axs[1].legend()

fig.tight_layout()

# Compare solved accelerations and numerically derived accelerations
fig2, ax2 = plt.subplots(1, 1)
Xdd_numerical = np.gradient(y[6, :], t)
plt.plot(tspan, Xdd[0,:], label=r"Solved $\ddot{x}_1$")
plt.plot(t, Xdd_numerical, label=r"Numerical $\ddot{x}_1$")
ax2.set_xlabel(r"Time ($s$)")
plt.title(r"Acceleration of Pendulum 1, $\ddot{x}_1$")
plt.grid(visible=True, which="major")
plt.legend()
fig2.tight_layout()

# Plotting the Lagrange Multiplicators
fig3, ax3 = plt.subplots(1, 1)
plt.plot(tspan, L[0, :], label=r"$\lambda_1$")
plt.plot(tspan, L[1, :], label=r"$\lambda_2$")
plt.plot(tspan, L[2, :], label=r"$\lambda_3$")
plt.plot(tspan, L[3, :], label=r"$\lambda_4$")
ax3.set_xlabel(r"Time ($s$)")
plt.title(r"Lagrangian Multiplicators/Factors")
plt.grid(visible=True, which="major")
plt.legend()
fig3.tight_layout()

plt.show()

sys.exit(0)
