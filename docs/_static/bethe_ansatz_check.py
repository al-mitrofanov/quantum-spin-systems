import numpy as np
from bethe_ansatz import bethe_one_mag, bethe_two_mag
from bethe_helpers import is_eigen, check_orthogonal


# load hamiltonians for magnons
number_of_spins = 10
exchange_constant = 1


hamiltonian_2_magnons = np.load('hamiltonian_2_magnons_10_spins.npy')

wavenumbers_one_magnon, states_one_magnon, energies_one_magnon = \
    bethe_one_mag(number_of_spins, exchange_constant)
wavenumbers_two_magnons, states_two_magnons, energies_two_magnons = \
    bethe_two_mag(number_of_spins, exchange_constant, hamiltonian_2_magnons)


# %% checking found solutions
for i in range(len(states_two_magnons)):
    if (is_eigen(hamiltonian_2_magnons, energies_two_magnons[i], states_two_magnons[i])):
        print("eigen")
    else:
        print("not eigen")

# from states_two_magnons to two_magnon_states
two_magnon_states = np.zeros((len(states_two_magnons[1]),len(states_two_magnons))).astype('complex')
for i in range(len(states_two_magnons)):
    two_magnon_states[:,i] = states_two_magnons[i].reshape(len(states_two_magnons[1]))

orthogonal = check_orthogonal(two_magnon_states)
