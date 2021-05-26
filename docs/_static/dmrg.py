import numpy as np
from hamiltonian_not_diagonal import form_not_diagonal_hamiltonian
from partial_trace_dmrg import partial_trace_initial, change_order

# %% 1 Constants
NUMBER_OF_SPINS_FINAL = 16
NUMBER_OF_SPINS_APPROXIMATION = 4
JX, JY, JZ = -1, -1, -1
SX = 0.5 * np.array([(0, 1), (1, 0)])
SY = 0.5 * np.array([(0, -1j), (1j, 0)])
SZ = 0.5 * np.array([(1, 0), (0, -1)])

# %% 2 create left and right Hamiltonian 2^m*2^m block
# 2.1 left
hamiltonian_left = -JZ * np.kron(SZ, SZ) - JX * np.kron(SX, SX) - JY * np.kron(SY, SY)
for i in range(NUMBER_OF_SPINS_APPROXIMATION - 2):
    i = i + 1
    sx_tilde_left = np.kron(np.eye(2 ** i), SX)
    sy_tilde_left = np.kron(np.eye(2 ** i), SY)
    sz_tilde_left = np.kron(np.eye(2 ** i), SZ)
    hamiltonian_left = (np.kron(hamiltonian_left, np.eye(2)) -
                        JZ * np.kron(sz_tilde_left, SZ) -
                        JX * np.kron(sx_tilde_left, SX) -
                        JY * np.kron(sy_tilde_left, SY)).astype(float)

# tilde operators to unite left and right blocks
sx_tilde_left = np.kron(np.eye(2), sx_tilde_left)
sy_tilde_left = np.kron(np.eye(2), sy_tilde_left)
sz_tilde_left = np.kron(np.eye(2), sz_tilde_left)

#  2.2 right
hamiltonian_right = -JZ * np.kron(SZ, SZ) - JX * np.kron(SX, SX) - JY * np.kron(SY, SY)
for i in range(NUMBER_OF_SPINS_APPROXIMATION - 2):
    i = i + 1
    sx_tilde_right = np.kron(SX, np.eye(2 ** i))
    sy_tilde_right = np.kron(SY, np.eye(2 ** i))
    sz_tilde_right = np.kron(SZ, np.eye(2 ** i))
    hamiltonian_right = (np.kron(np.eye(2), hamiltonian_right) -
                         JZ * np.kron(SZ, sz_tilde_right) -
                         JX * np.kron(SX, sx_tilde_right) -
                         JY * np.kron(SY, sy_tilde_right)).astype(float)

# tilde operators to unite left and right blocks
sx_tilde_right = np.kron(sx_tilde_right, np.eye(2))
sy_tilde_right = np.kron(sy_tilde_right, np.eye(2))
sz_tilde_right = np.kron(sz_tilde_right, np.eye(2))

# 2.3 Unite left and right Hamiltonians (superblock)
hamiltonian_superblock = (np.kron(hamiltonian_left, np.eye(2**NUMBER_OF_SPINS_APPROXIMATION)) +
                          np.kron(np.eye(2**NUMBER_OF_SPINS_APPROXIMATION), hamiltonian_right) -
                          JX * np.kron(sx_tilde_left, sx_tilde_right) -
                          JY * np.kron(sy_tilde_left, sy_tilde_right) -
                          JZ * np.kron(sz_tilde_left, sz_tilde_right)).astype(float)

# 2.4 Obtaining the ground state
eigen_vectors = np.linalg.eigh(hamiltonian_superblock)[1]

# 2.5 density matrixes of the left and right subsystems in the ground state
reduced_density_matrix_left, reduced_density_matrix_right = \
    partial_trace_initial(eigen_vectors[:, 0])
results = []
# %% 3 From new Hamiltonians in new basis
for i in range(NUMBER_OF_SPINS_APPROXIMATION + 1, NUMBER_OF_SPINS_FINAL + 1):
    # eigen vectors for each density matrix
    eigen_vectors_density_left = np.linalg.eigh(reduced_density_matrix_left, )[1]
    eigen_vectors_density_right = np.linalg.eigh(reduced_density_matrix_right)[1]
    # change order of eigen vectors
    eigen_vectors_density_left = change_order(eigen_vectors_density_left)
    eigen_vectors_density_right = change_order(eigen_vectors_density_right)
    # truncate eigen vectors
    eigen_vectors_density_left = eigen_vectors_density_left[:, :2**NUMBER_OF_SPINS_APPROXIMATION]
    eigen_vectors_density_right = eigen_vectors_density_right[:, :2**NUMBER_OF_SPINS_APPROXIMATION]
    # rewrite left and right hamiltonians in the new basis
    hamiltonian_left = eigen_vectors_density_left.T @  hamiltonian_left @ \
        eigen_vectors_density_left
    hamiltonian_right = eigen_vectors_density_right.T @  hamiltonian_right @ \
        eigen_vectors_density_right

    eigen_vectors_density_left = eigen_vectors_density_left[:2**NUMBER_OF_SPINS_APPROXIMATION,
                                                            :2**NUMBER_OF_SPINS_APPROXIMATION]
    eigen_vectors_density_right = eigen_vectors_density_right[:2**NUMBER_OF_SPINS_APPROXIMATION,
                                                              :2**NUMBER_OF_SPINS_APPROXIMATION]
    # rewrite tilde operators in eigen basis
    sx_tilde_left = np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION-1)), SX)
    sy_tilde_left = np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION-1)), SY)
    sz_tilde_left = np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION-1)), SZ)
    sx_tilde_right = np.kron(SX, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION-1)))
    sy_tilde_right = np.kron(SY, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION-1)))
    sz_tilde_right = np.kron(SZ, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION-1)))

    sx_tilde_left = eigen_vectors_density_left.T @  sx_tilde_left @ \
        eigen_vectors_density_left
    sy_tilde_left = eigen_vectors_density_left.T @  sy_tilde_left @ \
        eigen_vectors_density_left
    sz_tilde_left = eigen_vectors_density_left.T @  sz_tilde_left @ \
        eigen_vectors_density_left
    sx_tilde_right = eigen_vectors_density_right.T @  sx_tilde_right @ \
        eigen_vectors_density_right
    sy_tilde_right = eigen_vectors_density_right.T @  sy_tilde_right @ \
        eigen_vectors_density_right
    sz_tilde_right = eigen_vectors_density_right.T @  sz_tilde_right @ \
        eigen_vectors_density_right

    # generate left and right Hamiltonian with new site 2^m+1*2^m+1 each
    hamiltonian_left = (np.kron(hamiltonian_left, np.eye(2)) -
                    JZ * np.kron(sz_tilde_left, SZ) -
                    JX * np.kron(sx_tilde_left, SX) -
                    JY * np.kron(sy_tilde_left, SY)).astype(float)
    hamiltonian_right = (np.kron(np.eye(2), hamiltonian_right) -
                     JZ * np.kron(SZ, sz_tilde_right) -
                     JX * np.kron(SX, sx_tilde_right) -
                     JY * np.kron(SY, sy_tilde_right)).astype(float)
    # save eigen energies
    results.append([np.linalg.eigvals(hamiltonian_left)[0], np.linalg.eigvals(hamiltonian_left)[0]])
    # generate left and right tilde operators for superblock
    sx_tilde_left = np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION)), SX)
    sy_tilde_left = np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION)), SY)
    sz_tilde_left = np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION)), SZ)
    sx_tilde_right = np.kron(SX, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION)))
    sy_tilde_right = np.kron(SY, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION)))
    sz_tilde_right = np.kron(SZ, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION)))

    # generate superblock Hamiltonian 2^2m+2*2^2m+2
    hamiltonian_superblock = (np.kron(hamiltonian_left, np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION + 1))) +
                              np.kron(np.eye(2**(NUMBER_OF_SPINS_APPROXIMATION + 1)), hamiltonian_right) -
                              JX * np.kron(sx_tilde_left, sx_tilde_right) -
                              JY * np.kron(sy_tilde_left, sy_tilde_right) -
                              JZ * np.kron(sz_tilde_left, sz_tilde_right)).astype(float)

    # ground state of superblock
    eigen_vectors = np.linalg.eigh(hamiltonian_superblock)[1]
    reduced_density_matrix_left, reduced_density_matrix_right = \
        partial_trace_initial(eigen_vectors[:, 0])

print(results)
