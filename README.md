# CoM Dynamics (Double Compound Pendulum)
A simply exercise of modelling the motion of a double compound pendulum and computing the reaction forces in the base revolute joint.

The equation of motions were first modelled according to the formulation of the Lagrangian equations as shown below:

$$
\begin{bmatrix}M & {J_F}^T \\\ J_F & 0 \end{bmatrix}\begin{bmatrix}\ddot{x} \\\ \lambda\end{bmatrix} =
\begin{bmatrix}F \\\ -\dot{J}_F \dot{x} \end{bmatrix} \,
$$

where $M$ represents the mass matrix, $J_F$ the Jacobian of the constraint equations, and $F$ the external force vector. By integrating the equation of motions above, one obtains the positions and velocities in time (details in code).

Considering the classical example of a double compound pendulum, as a nice and simple validation step, we can generate the animation as shown below using the integrated positions.

<p align="center">
  <img src="https://github.com/RenZhen95/CoMMechanics/blob/main/animation/Lagrange1_alpha100beta100.gif">
</p>

A nice condition that we can use to check if the obtained solutions are valid is to ensure that the constraint equations are always zero (or at least close to zero).

<p align="center">
  <img src="https://github.com/RenZhen95/CoMMechanics/blob/main/plots/constraintEquations.png">
</p>

The reaction forces $R$ can be computed as $R = -{J_F}^T \lambda$, but alternatively we can also consider the CoM dynamics

$$
a_S = \frac{\sum m_i a_{S_i}}{\sum m_i}
$$

$$
\sum F^{ext} = \left( \sum m_i \right) \\ a_S ,
$$

which allows to obtain the reaction forces at the base revolute joint via

$$
\begin{bmatrix}F_x \\\ F_y\end{bmatrix} = \begin{bmatrix}M \ddot{x}_S \\\ M(\ddot{y}_S + g)\end{bmatrix}
$$

<p align="center">
  <img src="https://github.com/RenZhen95/CoMMechanics/blob/main/plots/reactionForces.png">
</p>
