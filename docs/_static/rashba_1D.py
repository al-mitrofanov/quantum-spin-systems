import numpy as np
import matplotlib.pyplot as plt

lx = 100 # length of the tight-binding sites
a = 1  # space constant
basis_kx = np.zeros(lx)
energy_1 = np.zeros(lx)
energy_2 = np.zeros(lx)
for i in range(lx):
    basis_kx[i] = -np.pi/a + (2*np.pi/(lx*a))*(i+1)
    energy_1[i] = (basis_kx[i]**2)/(2) - basis_kx[i]/4
    energy_2[i] = (basis_kx[i]**2)/(2) + basis_kx[i]/4

plt.plot(basis_kx, energy_1, label='Spin up')
plt.plot(basis_kx, energy_2, label='Spin down')
plt.legend()
plt.title('Dispersion law for 1D system in y basis')
plt.xlabel('k_x')
plt.ylabel('Energy')