import numpy as np

def add_spin_interaction(fist_spin, second_spin, new_hamiltonian, full_basis, j_fm_afm):
    list_double_count = []
    for i in range(len(full_basis)):
        if i not in list_double_count:
            if full_basis[i][fist_spin] == full_basis[i][second_spin]:
                new_hamiltonian[i, i] = new_hamiltonian[i, i] - j_fm_afm[2]/4
            else:
                new_hamiltonian[i, i] = new_hamiltonian[i, i] + j_fm_afm[2]/4
                looking_state = full_basis[i].copy()
                temp = looking_state[fist_spin]
                looking_state[fist_spin] = looking_state[second_spin]
                looking_state[second_spin] = temp
                l = find_state_in_the_basis(looking_state, full_basis)
                new_hamiltonian[i, l] = -0.5 * j_fm_afm[1]
                new_hamiltonian[l, i] = -0.5 * j_fm_afm[1]
                new_hamiltonian[l, l] = new_hamiltonian[l, l] + j_fm_afm[2]/4
                list_double_count.append(l)
    new_hamiltonian = new_hamiltonian.astype(np.float32)



def find_state_in_the_basis(looking_state, basis):
    for i in range(2 ** (len(looking_state))):
        if (basis[i] == looking_state).all():
            return i