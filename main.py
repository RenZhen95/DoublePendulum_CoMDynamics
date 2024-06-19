import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Geometry and inertia
m1 = 1   # kg
m2 = 1   # kg
g = 9.81 # m/s^2
l1 = 0.3 # m
l2 = 0.3 # m

# def getCOM():

# Initial conditions
r1_0 = np.array([np.cos(np.pi/3), -np.sin(np.pi/3), 0])
r1_0 = l1*r1_0

r2rel_0 = np.array([np.cos(np.pi/4), -np.sin(np.pi/4), 0])
r2rel_0 = l2*r2rel_0
r2_0 = r1_0 + r2rel_0

# # Quick plot to check initial conditions
# plt.plot([0, r1_0[0]], [0, r1_0[1]], 'o-')
# plt.plot([r1_0[0], r2_0[0]], [r1_0[1], r2_0[1]], 'o-')
# plt.show()

def getYd(t, Y_):
    X = Y_[0:5]
    Xd = Y_[5:10]


    Xdd_sub = np.matmul(
        np.matmul(
            J_T, np.linalg.inv(np.matmul(J, np.matmul(invM, J_T)))
        ),                                         # 5 x 3 matrix
        (np.matmul(J, np.matmul(invM, F)) + Jbg) # 3 x 1 matrix
    )
    Xdd = np.matmul(invM, F - Xdd_sub)

    return np.concatenate((Xd, Xdd))

# Integration
dt = 0.001
tspan = np.arange(0, 1.001, dt)
sol = solve_ivp(
    getYd, [0.0, 1.0], Y0, t_eval=tspan, rtol=1e-8, atol=1e-12
)

X_t = sol.y[0:5, :]
Xd_t = sol.y[5:, :]


