# CoM Dynamics (Double Compound Pendulum)
A simply exercise of modelling the motion of a double compound pendulum and computing the reaction forces in the base revolute joint.

The equation of motions were first modelled according to the formulation of the Lagrangian equations as shown below:

$$
\begin{bmatrix}M & {J_F}^T \\ J_F & 0 \end{bmatrix}\begin{bmatrix}\ddot{x} \\ \lambda\end{bmatrix} =
\begin{bmatrix}F \\ -\dot{J}_F \lambda \end{bmatrix}
$$