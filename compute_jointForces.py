import pickle
import os, sys
import numpy as np
import matplotlib.pyplot as plt

with open("integratedCoordinates.pkl", 'rb') as handle:
    y = pickle.load(handle)

# Geometry and inertia
m1 = 1           # kg
m2 = 1           # kg
g = 9.81         # m/s^2
l1 = 0.3         # m
l2 = 0.3         # m
J1 = m1*l1*l1/12 # kgm^2
J2 = m2*l2*l2/12 # kgm^2

# Integrated positions with time
X_t = y[0:6, :]
# Integrated velocities with time
Xd_t = y[6:12, :]

# CoM
M = m1+m2
# x-Coordinate
xS = (m1*X_t[0, :] + m2*X_t[3, :]) / M
# y-Coordinate
yS = (m1*X_t[1, :] + m2*X_t[4, :]) / M

# Get acceleration of CoM
xSdd = np.gradient(np.gradient(xS))
ySdd = np.gradient(np.gradient(yS))

# Reaction forces
Fx = M*xSdd
Fy = M*ySdd + M*g

# Plotting accelerations and joint force
title_text = "Parameters: "
title_text += r"$m_1:$ "
title_text += f"{m1} kg, "
title_text += r"$l_1:$ "
title_text += f"{l1} m | "
title_text += r"$m_2:$ "
title_text += f"{m2} kg, "
title_text += r"$l_2:$ "
title_text += f"{l2} m"
# x-axis
figX, axsX = plt.subplots(2, 1, sharex=True)
dt = 0.001
tspan = np.arange(0, 5.001, dt)
axsX[0].plot(tspan, xSdd)
axsX[0].set_ylabel(r"$\frac{m}{s^2}$")
axsX[0].set_title(r"Acceleration of CoM ($\ddot{x}_S$)", loc="left")
axsX[1].plot(tspan, Fx)
axsX[1].set_ylabel(r"$N$")
axsX[1].set_title(r"Reaction force ($F_x$)", loc="left")

axsX[0].grid(visible=True, which="major")
axsX[1].grid(visible=True, which="major")

axsX[1].set_xlabel(r"Time ($s$)")
figX.suptitle(title_text, x=0.115, y=0.95, ha="left")
figX.tight_layout()

# y-axis
figY, axsY = plt.subplots(2, 1, sharex=True)
dt = 0.001
tspan = np.arange(0, 5.001, dt)
axsY[0].plot(tspan, ySdd)
axsY[0].set_ylabel(r"$\frac{m}{s^2}$")
axsY[0].set_title(r"Acceleration of CoM ($\ddot{y}_S$)", loc="left")
axsY[1].plot(tspan, Fy)
axsY[1].set_ylabel(r"$N$")
axsY[1].set_title(r"Reaction force ($F_y$)", loc="left")

axsY[0].grid(visible=True, which="major")
axsY[1].grid(visible=True, which="major")

axsY[1].set_xlabel(r"Time ($s$)")
figY.suptitle(title_text, x=0.157, y=0.95, ha="left")
figY.tight_layout()
plt.show()

sys.exit(0)
