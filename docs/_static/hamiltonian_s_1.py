import numpy as np


def form_hamiltonian_s_1(exchange_constant, number_of_spins, borders_opened=True,
                             anisotropy=[1, 1, 1]):
    """
    Return Heisenberg Hamiltonian for n spins with S=1
    exchange_constant - exchange interaction constant, float
    number_of_spins - number of spins, int
    anisotropy - anisotropy constants, no anisotropy by defaults
    borders_opened - borders are opened if True, closed otherwise
    """
    jx = exchange_constant * anisotropy[0]
    jy = exchange_constant * anisotropy[1]
    jz = exchange_constant * anisotropy[2]
    SX_1 = 1/(np.sqrt(2)) * np.array([(0, 1, 0), (1, 0, 1), (0, 1, 0)])
    SY_1 = 1/(np.sqrt(2)) * np.array([(0, -1j, 0), (1j, 0, -1j), (0, 1j, 0)])
    SZ_1 = np.array([(1, 0, 0), (0, 0, 0), (0, 0, -1)])

    hamiltonian = -jz * np.kron(SZ_1, SZ_1) - jx * np.kron(SX_1, SX_1) - jy * np.kron(SY_1, SY_1)

    for i in range(number_of_spins - 2):
        i += 1
        sx_tilde = np.kron(np.eye(3 ** i), SX_1)
        sy_tilde = np.kron(np.eye(3 ** i), SY_1)
        sz_tilde = np.kron(np.eye(3 ** i), SZ_1)
        hamiltonian = np.kron(hamiltonian, np.eye(3)) - jz * np.kron(sz_tilde, SZ_1) - \
            jx * np.kron(sx_tilde, SX_1) - jy * np.kron(sy_tilde, SY_1)
    if borders_opened is False:
        sx_tilde = np.kron(SX_1, np.eye(3 ** (number_of_spins - 2)))
        sy_tilde = np.kron(SY_1, np.eye(3 ** (number_of_spins - 2)))
        sz_tilde = np.kron(SZ_1, np.eye(3 ** (number_of_spins - 2)))
        hamiltonian = hamiltonian - jz * np.kron(sz_tilde, SZ_1) - \
            jx * np.kron(sx_tilde, SX_1) - jy * np.kron(sy_tilde, SY_1)

    return hamiltonian.astype(float)



