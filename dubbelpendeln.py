import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

m = 1
l = 1
g = 9.81
theta1 = np.pi/10
theta2 = np.pi/10
p1 = 0
p2 = 0
y0 = np.array([theta1, theta2, p1, p2])

h = 0.1
t0 = 0
t1 = 10
tspan = (t0, t1)
tt = np.arange(t0, t1, h)


def dubbelpendel_ode(t, y):
    theta1, theta2, p1, p2 = y
    # variabler som underlättar uträkningen
    cos_diff = np.cos(theta1 - theta2) # täljaren
    cos_diff2 = np.cos(theta1 - theta1)**2  # nämnaren
    factor = 6/m*(l**2) 

    # vinkelhastigheter:
    theta1_deriv = factor * (((2*p1 - 3 * cos_diff)*p2) / (16 - 9*cos_diff2))
    theta2_deriv = factor * (((8*p2 - 3 * cos_diff)*p1) / (16 - 9*cos_diff2))

    # variabler som underlättar uträkningen
    sin_diff = np.sin(theta1 - theta2) 
    sin_theta1 = np.sin(theta1)
    sin_theta2 = np.sin(theta2)
    factor2 = -(1/2)*m*(l**2)

    # förändring av rörelsemängden:
    p1_deriv = factor2 * (((theta1_deriv * theta2_deriv) * sin_diff) + 3 * (g/l) * sin_theta1)
    p2_deriv = factor2 * (((-theta1_deriv * theta2_deriv) * sin_diff) + (g/l) * sin_theta2)

    return [theta1_deriv, theta2_deriv, p1_deriv, p2_deriv]


sol = solve_ivp(dubbelpendel_ode, tspan, y0, t_eval=tt)

plt.plot(sol.t, sol.y[0], label="y(t)")
plt.plot(sol.t, sol.y[1], label="y'(t)")
plt.title('Ett system av ODE:er, demo')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()