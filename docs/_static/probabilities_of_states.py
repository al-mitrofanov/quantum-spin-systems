import numpy as np

diag_basis = np.diag(np.ones(8))
ksi_1 = diag_basis[:, 0]
ksi_2 = diag_basis[:, 1]

psi = np.zeros([8, 1], dtype=complex)  # set the state
psi[0] = 0.3   # c_1 = 0.3
psi[1] = 0.9*1j   # c_2 = 0.9i
p_1 = ksi_1.conj().T @ psi @ psi.conj().T @ ksi_1  # probability of the 1st basis state
p_2 = ksi_2.conj().T @ psi @ psi.conj().T @ ksi_2  # probability of the 2nd basis state
print(p_1 == abs(psi[0])**2, p_2 == abs(psi[1])**2)  # compare with analytical results