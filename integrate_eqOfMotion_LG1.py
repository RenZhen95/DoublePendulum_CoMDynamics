import sys
import pickle
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.integrate import solve_ivp

# Geometry and inertia
m1 = 1           # kg
m2 = 1           # kg
g = 9.81         # m/s^2
l1 = 0.3         # m
l2 = 0.3         # m
J1 = m1*l1*l1/12 # kgm^2
J2 = m2*l2*l2/12 # kgm^2

# Mass matrix
M = np.diag([m1, m1, J1, m2, m2, J2])

# Force array
F = np.array([0, -m1*g, 0, 0, -m2*g, 0])

# Initial conditions
# - Position
X0 = np.zeros(6)
X0[2] = 85 * (np.pi/180) # 85
X0[0] = (l1/2) * np.sin(X0[2])
X0[1] = -(l1/2) * np.cos(X0[2])

X0[5] = 70 * (np.pi/180) # 70
X0[3] = (l1*np.sin(X0[2])) + (l2/2*np.sin(X0[5]))
X0[4] = -(l1*np.cos(X0[2])) - (l2/2*np.cos(X0[5]))

# - Velocity
Xd0 = np.zeros(6)

# Initializing the state vector
Y0 = np.concatenate((X0, Xd0))

def getJacobian(x):
    x1 = x[0]; y1 = x[1]; phi1 = x[2]
    x2 = x[3]; y2 = x[4]; phi2 = x[5]

    s1 = np.sin(phi1); c1 = np.cos(phi1)
    s2 = np.sin(phi2); c2 = np.cos(phi2)

    J = np.array([
        [2*x1, 2*y1, 0, 0, 0, 0],
        [c1, s1, -x1*s1 + y1*c1, 0, 0, 0],
        [0, 0, -l1*c1, 1, 0, -l2/2*c2],
        [0, 0, -l1*s1, 0, 1, -l2/2*s2]
    ])
    return J
def getJacobianDot(x, xd):
    x1 = x[0]; y1 = x[1]; phi1 = x[2]
    x2 = x[3]; y2 = x[4]; phi2 = x[5]

    x1d = xd[0]; y1d = xd[1]; phi1d = xd[2]
    x2d = xd[3]; y2d = xd[4]; phi2d = xd[5]

    s1 = np.sin(phi1); c1 = np.cos(phi1)
    s2 = np.sin(phi2); c2 = np.cos(phi2)

    Jd = np.array([
        [2*x1d, 2*y1d, 0, 0, 0, 0],
        [-phi1d*s1, phi1d*c1, -x1*phi1d*c1 - x1d*s1 - y1*phi1d*s1 + y1d*c1, 0, 0, 0],
        [0, 0, l1*phi1d*s1, 0, 0, l2/2*phi2d*s2],
        [0, 0, -l1*phi1d*c1, 0, 0, -l2/2*phi2d*c2]
    ])
    return Jd
def getConstraints(x):
    x1 = x[0]; y1 = x[1]; phi1 = x[2]
    x2 = x[3]; y2 = x[4]; phi2 = x[5]

    s1 = np.sin(phi1); c1 = np.cos(phi1)
    s2 = np.sin(phi2); c2 = np.cos(phi2)

    f = np.array([
        x1*x1 + y1*y1 - (l1*l1)/4,
        x1*c1 + y1*s1,
        x2 - l1*s1 - l2/2*s2,
        y2 + l1*c1 + l2/2*c2
    ])
    return f
def getConstraintsDot(x, xd):
    x1 = x[0]; y1 = x[1]; phi1 = x[2]
    x2 = x[3]; y2 = x[4]; phi2 = x[5]

    x1d = xd[0]; y1d = xd[1]; phi1d = xd[2]
    x2d = xd[3]; y2d = xd[4]; phi2d = xd[5]

    s1 = np.sin(phi1); c1 = np.cos(phi1)
    s2 = np.sin(phi2); c2 = np.cos(phi2)

    fdot = np.array([
        2*x1*x1d + 2*y1*y1d,
        -x1*phi1d*s1 + x1d*c1 + y1*phi1d*c1 + y1d*s1,
        x2d - l1*phi1d*c1 - l2/2*phi2d*c2,
        y2d - l1*phi1d*s1 - l2/2*phi2d*s2
    ])
    return fdot

# Baumgarte-Stabilization Parameters
alpha = 100
beta = 100

# Integration
dt = 0.001
tspan = np.arange(0, 5.001, dt)

# Reaction forces
R = []
# Lambdas (Lagrange multiplicators)
L = []
# Accelerations
Xddarray = []

def getYd(t, Y_):
    X = Y_[0:6]
    Xd = Y_[6:12]

    J = getJacobian(X)
    Jd = getJacobianDot(X, Xd)
    f = getConstraints(X)
    fd = getConstraintsDot(X, Xd)

    Jbg = np.matmul(Jd, Xd) + alpha*fd + beta*f
    J_T = np.transpose(J)
    invM = np.linalg.inv(M)

    # Lagrange multiplicator
    Lm = np.matmul(
        np.linalg.inv(np.matmul(J, np.matmul(invM, J_T))),
        (np.matmul(J, np.matmul(invM, F)) + Jbg)
    )
    L.append(Lm)

    # Reaction forces
    R.append(np.matmul(J_T, Lm))

    Xdd = np.matmul(invM, F - np.matmul(J_T, Lm))
    Xddarray.append(Xdd)

    return np.concatenate((Xd, Xdd))

# Integration
sol = solve_ivp(
    getYd, [0.0, 5.0], Y0, t_eval=tspan, rtol=1e-8, atol=1e-12
)

X_t = sol.y[0:6, :]
Xd_t = sol.y[6:, :]

to_save = {
    "t_eval": sol.t,
    "integratedCoordinates": sol.y,
    "lambdas": L,
    "reactionForces": R,
    "Xdd": Xddarray
}
with open("savedResults.pkl", 'wb') as handle:
    pickle.dump(to_save, handle)

# Check that constraint equations are always zero
F_t = getConstraints(X_t)
fig1, axs1 = plt.subplots(1, 1)
axs1.plot(tspan, F_t[0,:], label=r"$f_1 (t)$")
axs1.plot(tspan, F_t[1,:], label=r"$f_2 (t)$")
axs1.plot(tspan, F_t[2,:], label=r"$f_3 (t)$")
axs1.plot(tspan, F_t[3,:], label=r"$f_4 (t)$")
axs1.grid(visible=True, which="major")
axs1.set_ylabel(r"$f$", rotation=0, y=0.9, labelpad=9.0)
axs1.set_xlabel(r"Time ($s$)")
axs1.set_title("Progression of constraint equations", loc="left")
axs1.legend()
fig1.tight_layout()

# Animation
fps = 175
# [Time per Frame / Timestep] = [Steps/Frame]
step_perFrame = int(1/fps/dt)

def make_plot(t_i, _frames):
    # CoM of Pendulum 1
    x1 = X_t[0, t_i]
    y1 = X_t[1, t_i]
    # CoM of Pendulum 2
    x2 = X_t[3, t_i]
    y2 = X_t[4, t_i]

    # End of Pendulum 2
    end2x = x1*2 + (l2*np.sin(X_t[5, t_i]))
    end2y = y1*2 - (l2*np.cos(X_t[5, t_i]))

    # Plot
    ax.plot([0, x1*2], [0, y1*2], 'o-', markersize=12)
    ax.plot([x1*2, end2x], [y1*2, end2y], 'o-', markersize=12)
    ax.plot([x1], [y1], 'o', markersize=12, color='black')
    ax.text(x1+0.02, y1, r'$m_1$')
    ax.plot([x2], [y2], 'o', markersize=12, color='black')
    ax.text(x2+0.02, y2, r'$m_2$')
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlim([-2.05*l1, 2.05*l1])
    ax.set_ylim([-2.05*l1, 0.25*l1])
    plt.grid(visible=True)
    plt.savefig('frames/_img{:04d}.png'.format(t_i//step_perFrame), dpi=72)
    image = imageio.imread('frames/_img{:04d}.png'.format(t_i//step_perFrame))
    _frames.append(image)
    plt.cla()

fig, ax = plt.subplots(1,1, figsize=(9,9))

frames = []
for i in range(0, int((tspan.size)/2), step_perFrame):
    make_plot(i, frames)

imageio.mimsave(f"Lagrange1_alpha{alpha}beta{beta}.gif", frames, loop=1)

fig1.show()
input("Press ENTER to quit ... ")

sys.exit(0)
