import numpy as np

x_eigen = np.matrix([[0.707, 0.707], [0.707, -0.707]])
y_eigen = np.matrix([[0.707, 0.707], [1j*0.707, -1j*0.707]])
z_eigen = np.matrix([[1, 0], [0, 1]])
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

sigma_y_in_x_basis = x_eigen.H @ sigma_y @ x_eigen
sigma_y_in_z_basis = z_eigen.H @ sigma_y @ z_eigen

