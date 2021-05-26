import numpy as np

def partial_trace_a(state):
    ro = state @ state.T
    ksi = np.eye(4)  # set of basis states
    ro_reduced = np.zeros([2])
    for i in range(4):
        ksi_full = np.kron(np.eye(2), np.matrix(ksi[:, i]).T)
        ro_reduced = ro_reduced + ksi_full.T @ ro @ ksi_full
    return ro_reduced


def partial_trace_b(state):
    ro = state @ state.T
    ksi = np.eye(2)
    ksi_full = np.zeros([8, 2, 4])
    ro_reduced = np.zeros([2])
    ksi_full = np.zeros([8, 2, 4])
    ksi_full[:, :, 0] = np.kron(np.array([[1], [0]]), np.kron(np.eye(2), np.array([[1], [0]])))
    ksi_full[:, :, 1] = np.kron(np.array([[1], [0]]), np.kron(np.eye(2), np.array([[0], [1]])))
    ksi_full[:, :, 2] = np.kron(np.array([[0], [1]]), np.kron(np.eye(2), np.array([[1], [0]])))
    ksi_full[:, :, 3] = np.kron(np.array([[0], [1]]), np.kron(np.eye(2), np.array([[0], [1]])))
    for i in range(4):
        ro_reduced = ro_reduced + ksi_full[:, :, i].T @ ro @ ksi_full[:, :, i]
    return ro_reduced


def partial_trace_c(state):
    ro = state @ state.T
    ksi = np.eye(4)  # set of basis states
    ro_reduced = np.zeros([2])
    for i in range(4):
        ksi_full = np.kron(np.matrix(ksi[:, i]).T, np.eye(2))
        ro_reduced = ro_reduced + ksi_full.T @ ro @ ksi_full
    return ro_reduced