import numpy as np

def partial_trace_initial(state):
    if state.shape == (len(state),):
        ro = np.matrix(state).T @ np.matrix(state)  # density matrix for the state
    else:
        ro = state @ state.T
    m = int(np.log2(len(state))/2)  # number of spins
    ksi = np.eye(2**m)  # set of basis vectors
    ro_reduced_left = np.zeros([2**m, 2**m])
    ro_reduced_right = np.zeros([2**m, 2**m])
    for i in range(2**m):
        ksi_full_left = np.kron(np.eye(2**m), np.matrix(ksi[:, i]).T)
        ksi_full_right = np.kron(np.matrix(ksi[:, i]).T, np.eye(2**m))
        ro_reduced_left = ro_reduced_left + ksi_full_left.T @ ro @ ksi_full_left
        ro_reduced_right = ro_reduced_right + ksi_full_right.T @ ro @ ksi_full_right
    return ro_reduced_left, ro_reduced_right


def change_order(eigen_vectors):
    eigen_vectors_new = np.zeros([len(eigen_vectors), len(eigen_vectors)])
    for i in range(len(eigen_vectors)):
        eigen_vectors_new[:, i] = (eigen_vectors[:, -i-1]).reshape(len(eigen_vectors))
    return eigen_vectors_new
