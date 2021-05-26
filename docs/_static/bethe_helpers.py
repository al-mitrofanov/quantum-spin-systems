import numpy as np
import math
from scipy.optimize import fsolve
from scipy.optimize import newton_krylov
from scipy.optimize import anderson


def clean(A, order):
    # A =np.matrix(A)
    # A= complex(A)
    if len(A.shape) > 1:
        length = A.shape[0]
        width = A.shape[1]
        for j in range(width):
            for i in range(length):
                if abs(A[i, j].real) < order:
                    A[i, j] = 1j * A[i, j].imag
                if abs(A[i, j].imag) < order:
                    A[i, j] = A[i, j].real
        return A
    else:
        length = A.shape[0]
        for i in range(length):
            if abs(A[i].real) < order:
                A[i] = 1j * A[i].imag
            if abs(A[i].imag) < order:
                A[i] = A[i].real
        return A


def is_eigen(H, E, psi):
    # IfEigen(H_diag, eig_val[0], eig_vec[:, 0])
    PsiE1 = H.dot(psi)  # matrix multiplication of Matrix H and vector psi
    PsiE2 = E * psi  # normalyzed eigen states
    PsiE1 = np.round(PsiE1, 3)
    PsiE2 = np.round(PsiE2, 3)
    sum_elem = 0  # sum of elements
    for i in range(PsiE1.shape[0]):
        if PsiE1[i] == PsiE2[i]:
            sum_elem = sum_elem + 1
    if sum_elem == PsiE1.shape[0]:
        return True
    else:
        return False


def normalize(a):
    norm = 0
    if len(a.shape) > 1:
        width = a.shape[0]
        length = a.shape[1]
        a_normalyzed = np.zeros([width, length])
        a_normalyzed = a_normalyzed.astype(complex)
        for j in range(length):
            for i in range(width):
                norm = norm + abs(a[i, j]) ** 2
            a_normalyzed[:, j] = (1 / math.sqrt(norm)) * a[:, j]
            norm = 0
        return a_normalyzed
    else:
        width = a.shape[0]
        a_normalyzed = np.zeros([width, 1])
        a_normalyzed = a_normalyzed.astype(complex)
        for i in range(width):
            norm = norm + abs(a[i]) ** 2
        a_normalyzed = (1 / math.sqrt(norm)) * a
        return a_normalyzed


def delta_function(a, b):
    if a == b:
        result = 1
    else:
        result = 0
    return result


def change_order_matrix(matrix_initial):
    matrix_changed_order = np.zeros([matrix_initial.shape[0],
                                     matrix_initial.shape[1]]).astype(complex)
    for i in range(matrix_initial.shape[0]):
        matrix_changed_order[i, :] = matrix_initial[matrix_initial.shape[0] - i - 1, :]
    return matrix_changed_order


def check_orthogonal(vectorList):
    size = len(vectorList)
    res = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            res[i, j] = np.real(np.matmul(np.transpose(np.conjugate(vectorList[i])), vectorList[j])) ** 2
            if (res[i][j] < 1e-15):
                res[i][j] = 0
    return res


def orthogonal(vec1, vec2):
    if (np.real(np.matmul(np.transpose(np.conjugate(vec1)), vec2)) ** 2) < 1e-10:
        return True
    else:
        return False
