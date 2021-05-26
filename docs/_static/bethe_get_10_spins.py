import numpy as np
from hamiltonian_not_diagonal import form_not_diagonal_hamiltonian
from hamiltonian_diagonal import form_heisenberg_block_diag

# load hamiltonians for magnons
number_of_spins = 10
exchange_constant = 1

# to save block of hamiltonian with 2 magnons
hamiltonian = form_not_diagonal_hamiltonian(exchange_constant, number_of_spins, borders_opened=False)
hamiltonian_diagonal = form_heisenberg_block_diag(hamiltonian)
hamiltonian_2_magnons = hamiltonian_diagonal[0][11:56, 11:56]
np.save('hamiltonian_2_magnons_10_spins2', hamiltonian_2_magnons)

