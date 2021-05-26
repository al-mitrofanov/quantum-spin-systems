import numpy as np
from hamiltonian_not_diagonal import *
from add_interaction import *


number_of_spins = 2
two_dimensional_coupling = True
couple_two_lines = False
# exchange anisotropy of the 1st_ham, 2nd and between
j = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])

hamiltonian_1 = form_not_diagonal_hamiltonian(1, number_of_spins)[0]
hamiltonian_2 = form_not_diagonal_hamiltonian(1, number_of_spins)[0]

full_basis = get_basis_vectors(2 * number_of_spins)
new_hamiltonian = np.kron(hamiltonian_1, np.identity(2**number_of_spins)) + \
    np.kron(np.identity(2**number_of_spins), hamiltonian_2)

if two_dimensional_coupling:
    for i in range(number_of_spins):
        first_spin = i
        second_spin = i + number_of_spins
        add_spin_interaction(first_spin, second_spin, new_hamiltonian, full_basis, -j[2])


if couple_two_lines:
    add_spin_interaction(number_of_spins-1, number_of_spins, new_hamiltonian, full_basis, j[2])
    hamiltonian_3 = form_not_diagonal_hamiltonian(1, 2*number_of_spins)[0]

