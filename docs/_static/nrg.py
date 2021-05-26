import numpy as np
from hamiltonian_not_diagonal import form_not_diagonal_hamiltonian

NUMBER_OF_SPINS_FINAL = 12
NUMBER_OF_SPINS_APPROXIMATION = 4
JX, JY, JZ = 1, 1, 1
SX = 0.5 * np.array([(0, 1), (1, 0)])
SY = 0.5 * np.array([(0, -1j), (1j, 0)])
SZ = 0.5 * np.array([(1, 0), (0, -1)])

# %% Create Hamiltonian m*m by the iterative way

hamiltonian = -JZ * np.kron(SZ, SZ) - JX * np.kron(SX, SX) - JY * np.kron(SY, SY)
for i in range(NUMBER_OF_SPINS_APPROXIMATION - 2):
    i = i + 1
    sx_tilde = np.kron(np.eye(2 ** i), SX)
    sy_tilde = np.kron(np.eye(2 ** i), SY)
    sz_tilde = np.kron(np.eye(2 ** i), SZ)
    hamiltonian = np.kron(hamiltonian, np.eye(2)) - JZ * np.kron(sz_tilde, SZ) - \
        JX * np.kron(sx_tilde, SX) - JY * np.kron(sy_tilde, SY)
# % tilde operators
sx_tilde = np.kron(np.eye(2), sx_tilde)
sy_tilde = np.kron(np.eye(2), sy_tilde)
sz_tilde = np.kron(np.eye(2), sz_tilde)

# %% Normal Renormalization Group

for i in range(NUMBER_OF_SPINS_APPROXIMATION+2, NUMBER_OF_SPINS_FINAL+2):
    [eigen_values, eigen_vectors] = np.linalg.eigh(hamiltonian)
    hamiltonian = eigen_vectors.T @ hamiltonian @ eigen_vectors
    hamiltonian = hamiltonian[:2**NUMBER_OF_SPINS_APPROXIMATION,
                              :2**NUMBER_OF_SPINS_APPROXIMATION]
    eigen_vectors = eigen_vectors[:2**NUMBER_OF_SPINS_APPROXIMATION,
                                  :2**NUMBER_OF_SPINS_APPROXIMATION]
    sx_tilde = eigen_vectors.T @ sx_tilde @ eigen_vectors
    sy_tilde = eigen_vectors.T @ sy_tilde @ eigen_vectors
    sz_tilde = eigen_vectors.T @ sz_tilde @ eigen_vectors
    hamiltonian = np.kron(hamiltonian, np.eye(2)) - JZ * np.kron(sz_tilde, SZ) - \
        JX * np.kron(sx_tilde, SX) - JY * np.kron(sy_tilde, SY)

[eigen_values_nrg_hamiltonian, eigen_vectors_nrg_hamiltonian] = np.linalg.eigh(hamiltonian)

# %% Compare with certain results

[eigen_values_nrg_certain, eigen_vectors_nrg_certain] = \
    np.linalg.eigh(form_not_diagonal_hamiltonian(1, NUMBER_OF_SPINS_FINAL)[0])
