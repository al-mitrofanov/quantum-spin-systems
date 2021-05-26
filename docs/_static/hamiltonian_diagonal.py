import numpy as np

def form_heisenberg_block_diag(hamiltonian_not_diag):
    """Return Heisenberg Hamiltonian in the block-diagonal form"""
    basis_list_not_diagonal = hamiltonian_not_diag[1]
    hamiltonian_matrix = hamiltonian_not_diag[0]
    number_of_spins = int(np.log2(len(hamiltonian_not_diag[0])))

    basis_list_diag = get_basis_diag(number_of_spins)  # getting diagonal basis
    eye_like_matrix = np.zeros((2 ** number_of_spins, 2 ** number_of_spins)) # form matrix
    filled_elements = list(np.arange(len(hamiltonian_matrix)))
    for i in range(2 ** number_of_spins):
        for j in filled_elements:
            bool_array = basis_list_not_diagonal[i] == basis_list_diag[j]
            if np.all(bool_array):
                eye_like_matrix[i, j] = 1
                filled_elements.remove(j)
                break
    # change basis order
    hamiltonian_matrix = eye_like_matrix.T @ hamiltonian_matrix @ eye_like_matrix

    return hamiltonian_matrix, basis_list_diag


def get_basis_diag(number_of_spins):
    counter_of_states = 0
    basis_diagonal = []
    basis_diagonal.append(np.array([1]*number_of_spins))

    for block in range(number_of_spins):
        block += 1
        number_of_states_in_block = int(np.math.factorial(number_of_spins) /
                                        (np.math.factorial(number_of_spins - block) *
                                          np.math.factorial(block)))
        block_basis = []
        for j in range(number_of_states_in_block):
            counter_of_states += 1
            second_counter = 0
            while True:
                # get rid of '0b'
                k_initial = bin(2 ** number_of_spins - 2 ** (block) - second_counter)[2:]
                # set size to the number of spins
                k_initial = '0'*(number_of_spins-len(k_initial)) + k_initial
                state = np.array([int(i) for i in k_initial])  # fron str to int
                number_of_zeros = number_of_spins - np.count_nonzero(state)
                if not was_it_before(state, block_basis) and number_of_zeros == block:
                    block_basis.append(state)
                    break
                second_counter += 1
        for i in block_basis:
            basis_diagonal.append(i)
    return basis_diagonal

def was_it_before(state, block_basis):
    """ is state in block_basis"""
    for i in block_basis:
        bool_result = state == i
        if np.all(bool_result) == True:
            return True
    return False