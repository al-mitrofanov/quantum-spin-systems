import numpy as np


def square_operator_expected_value(state):
    """Return eigen value of square operator for 2 spin state"""
    sx_square = (np.kron(SX, np.eye(2)) + (np.kron(np.eye(2), SX)))
    sy_square = (np.kron(SY, np.eye(2)) + (np.kron(np.eye(2), SY)))
    sz_square = (np.kron(SZ, np.eye(2)) + (np.kron(np.eye(2), SZ)))
    square_operator = sx_square @ sx_square + sy_square @ sy_square + sz_square @ sz_square
    return state.T @ square_operator @ state


SX = 0.5 * np.array([(0, 1), (1, 0)])
SY = 0.5 * np.array([(0, -1j), (1j, 0)])
SZ = 0.5 * np.array([(1, 0), (0, -1)])

hamiltonian = np.kron(SX, SX) + np.kron(SY, SY) + np.kron(SZ, SZ)  # forming hamiltonian for 2 spins

[eigen_values, eigen_vectors] = np.linalg.eig(hamiltonian)  # eigen values and eigen vectors

print(*[square_operator_expected_value(eigen_vectors[:, i]) for i in range(4)])