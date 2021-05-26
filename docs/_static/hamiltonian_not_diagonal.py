import numpy as np


def form_not_diagonal_hamiltonian(exchange_constant, number_of_spins, borders_opened=True,
                                  anisotropy=[1, 1, 1]):
    """
    Return Heisenberg Hamiltonian for n spins
    exchange_constant - exchange interaction constant, float
    number_of_spins - number of spins, int
    anisotropy - anisotropy constants, no anisotropy by defaults
    borders_opened - borders are opened if True, closed otherwise
    """
    jx = exchange_constant * anisotropy[0]
    jy = exchange_constant * anisotropy[1]
    jz = exchange_constant * anisotropy[2]
    SX = 0.5 * np.array([(0, 1), (1, 0)])
    SY = 0.5 * np.array([(0, -1j), (1j, 0)])
    SZ = 0.5 * np.array([(1, 0), (0, -1)])

    hamiltonian = -jz * np.kron(SZ, SZ) - jx * np.kron(SX, SX) - jy * np.kron(SY, SY)

    for i in range(number_of_spins - 2):
        i += 1
        sx_tilde = np.kron(np.eye(2 ** i), SX)
        sy_tilde = np.kron(np.eye(2 ** i), SY)
        sz_tilde = np.kron(np.eye(2 ** i), SZ)
        hamiltonian = np.kron(hamiltonian, np.eye(2)) - jz * np.kron(sz_tilde, SZ) - \
            jx * np.kron(sx_tilde, SX) - jy * np.kron(sy_tilde, SY)
    if borders_opened is False:
        sx_tilde = np.kron(SX, np.eye(2 ** (number_of_spins - 2)))
        sy_tilde = np.kron(SY, np.eye(2 ** (number_of_spins - 2)))
        sz_tilde = np.kron(SZ, np.eye(2 ** (number_of_spins - 2)))
        hamiltonian = hamiltonian - jz * np.kron(sz_tilde, SZ) - \
            jx * np.kron(sx_tilde, SX) - jy * np.kron(sy_tilde, SY)

    return hamiltonian.astype(float), get_basis_vectors(number_of_spins)

def get_basis_vectors(number_of_spins):
    """ Return list of n basis states, in the tensor_product order """
    basis_list = []
    for i in range(2 ** number_of_spins):
        k = bin(2 ** number_of_spins - 1 - i)
        vector = np.zeros(number_of_spins).astype(int)
        for i in range(len(k) - 2):
            vector[-i - 1] = int(k[-i - 1])
        basis_list.append(vector)

    return basis_list
