import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

hbar = 1.05e-34
q = 1.6e-19

hamiltonian = q * np.array([[1, 0], [0, -1]])  # gamma*B*hbar = q
omega = 2 * q / hbar  # Rabi frequency
sx_state = 1/np.sqrt(2)*np.array([[1], [1]])
s_minus_x_state = 1/np.sqrt(2)*np.array([[1], [-1]])

psi0 = sx_state  # initial state along the x axis
time = 4  # fs
time_step = 0.05  # fs
t = 0

psi_t = []  # create a variable
expect_x = np.zeros((int(time / time_step)), dtype=float)  # create a variable
expect_minus_x = np.zeros((int(time / time_step)), dtype=float)  # create a variable

for j in range(int(time / time_step)):
    U = expm(-1j * hamiltonian * t / hbar)  # evolution operator
    psi_t.append(U @ psi0)  # Calculate psi for every time step
    # expectation to get x polarization after a measurement
    expect_x[j] = (abs(sx_state.T @ psi_t[j])) ** 2
    # expectation to get -x polarization after a measurement
    expect_minus_x[j] = (abs(s_minus_x_state.T @ psi_t[j])) ** 2
    t = t + time_step * 1e-15

plt.figure(1)
t = 1e-15 * np.linspace(0, time, int(time / time_step))
cos2 = np.cos(omega * t / 2) ** 2  # analytical solution
plt.plot(t, expect_x, label='expectation to get x')
plt.plot(t, expect_minus_x, label='expectation to get -x')
plt.plot(t, cos2, 'r^', label='analytical solution')
plt.xlabel('time, fs')
plt.ylabel('Probability')
plt.legend()
