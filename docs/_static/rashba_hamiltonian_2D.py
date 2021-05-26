import numpy as np
import matplotlib.pyplot as plt
import math
from bethe_ansatz import clean

# %% Rashba Hamiltonian

alpha = 1.5   # constant of RSO Hamiltonian
lx = 20  # length of the tight-binding sites
ly = 20  # width of the tight-binding sites
a = 1  # lattice constant
m = 1  # effective mass
spin_basis = 'z'  # can be z or y
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=int)

basis_kx = np.arange(-np.pi / a, np.pi / a, 2 * np.pi / (lx * a), dtype=float)
basis_ky = np.arange(-np.pi / a, np.pi / a, 2 * np.pi / (ly * a), dtype=float)

hamRSOI = np.zeros([lx * ly * 2, lx * ly * 2], dtype=complex)
count = 0
if spin_basis == 'y':  # z basis
    for i in range(lx):
        for j in range(ly):
            hamRSOI[count, count + 1] = alpha * (-1j * basis_kx[i] + basis_ky[j])
            hamRSOI[count + 1, count] = alpha * (1j * basis_kx[i] + basis_ky[j])
            hamRSOI[count, count] = (basis_kx[i] ** 2 + basis_ky[j] ** 2) / (2 * m)
            hamRSOI[count + 1, count + 1] = (basis_kx[i] ** 2 + basis_ky[j] ** 2) / (2 * m)
            count += 2
else:  # y basis
    for i in range(lx):
        for j in range(ly):
            hamRSOI[count, count + 1] = alpha * (basis_ky[j])
            hamRSOI[count + 1, count] = alpha * (basis_ky[j])
            hamRSOI[count, count] = (basis_kx[i] ** 2 + basis_ky[j] ** 2) / (2 * m) + \
                alpha * (basis_kx[i])
            hamRSOI[count + 1, count + 1] = (basis_kx[i] ** 2 + basis_ky[j] ** 2) / (2 * m) - \
                alpha * (basis_kx[i])
            count += 2

eig_val, eig_vec = np.linalg.eig(hamRSOI)
eig_val = clean(eig_val, 1e-10)  # delete zeo values for clarity

# %% Properties of the every state, saved in band_split

count = 0
band_split = np.zeros([lx * ly, 2], dtype=float)
for i in range(lx):
    for j in range(ly):
        band_split[count, 0] = basis_kx[i]
        band_split[count, 1] = basis_ky[j]
        count += 1
band_split = np.kron(band_split, np.array([[1], [1]]))
band_split = np.hstack((band_split, np.zeros([2 * lx * ly, 4])))

spin_state = np.zeros([2 * lx * ly, 2], dtype=complex)
count = 0
for state_number in range(lx * ly * 2):
    eig_vec_line = eig_vec[:, state_number]
    if state_number % 2 == 0:
        count = 0
    else:
        count = 1
    for i in range(lx * ly * 2):
        if abs(eig_vec_line[i]) != 0:
            if abs(eig_vec_line[i]) != 1:
                band_split[i + count, 2] = eig_val[state_number]
                spin_state[i + count, 0] = eig_vec_line[i]
                spin_state[i + count, 1] = eig_vec_line[i + 1]

            else:
                band_split[i, 2] = eig_val[state_number]
                if count == 0:
                    spin_state[i, 0] = eig_vec_line[i]
                    spin_state[i, 1] = 0
                else:
                    spin_state[i, 0] = 0
                    spin_state[i, 1] = eig_vec_line[i]
            break

for i in range(2 * lx * ly):
    # add expectation values of spin operators
    psi = np.matrix([[spin_state[i, 0]], [spin_state[i, 1]]])
    ro = np.matmul(psi, psi.H)
    band_split[i, 3] = np.trace(np.matmul(0.5 * sigma_x, ro))
    band_split[i, 4] = np.trace(np.matmul(0.5 * sigma_y, ro))
    band_split[i, 5] = np.trace(np.matmul(0.5 * sigma_z, ro))

# separating bands, if needed
band1 = np.zeros([lx * ly, 6], dtype=float)
band2 = np.zeros([lx * ly, 6], dtype=float)
count = 0
for i in range(2 * lx * ly):
    if i % 2 == 0:
        band1[count, :] = band_split[i, 0:6]
    else:
        band2[count, :] = band_split[i, 0:6]
        count += 1

