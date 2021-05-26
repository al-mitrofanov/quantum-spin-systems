import numpy as np
from bethe_helpers import normalize, change_order_matrix, clean, is_eigen, orthogonal, delta_function
from scipy.optimize import fsolve
from scipy.optimize import newton_krylov
from scipy.optimize import anderson


def bethe_one_mag(n, J):
    """ Return eigen vectors with 1 magnon, wavenumbers and eigen energies
    n - number of spins
    J - exhcange constant
    """
    k = np.arange(0, 2 * np.pi, 2 * np.pi/n)  # k is a wavenumber
    states, eigen_energies = [], []
    ground_energy = -J * n/4  # ground state energy for FM coupling
    for i in range(n):
        state = np.zeros((n, 1)).astype('complex')
        for j in range(n):
            state[j, 0] = np.exp(-1j * k[i] * j)  # every state is a single wave with wavenumber k
        states.append(normalize(clean(state, 1e-10)))
        eigen_energies.append((1 - np.cos(k[i])) + ground_energy)  #
    return k, states, eigen_energies


def bethe_two_mag(n, J, hamiltonian_2magnons):
    """ Return eigen vectors with 2 magnons, 2 wavenumbers and eigen energies
    n - number of spins
    J - exhcange constant
    H - hamiltonian of the block with 2 magnons
    """
    # class C1: k1 is variable, k2 = 0 is constant
    k, states_result, energies = [], [], []
    ground_energy = -J * n/4
    length = int(n * (n - 1) / 2)  # amount of 2 magnon states

    for i in range(n):
        k2 = 2 * np.pi * i / n
        k.append(np.array([0, k2]))
        energies.append(J * (1 - np.cos(k2)) + ground_energy)
        count = 0
        state = np.zeros((length, 1)).astype('complex')
        for n1 in range(1, n):
            for n2 in range(n1 + 1, n + 1):
                state[count, 0] = np.exp(1j * (k2 * n2)) + np.exp(1j * (k2 * n1))
                count += 1
        states_result.append(change_order_matrix(clean(normalize(state), 1e-10)))

    # class C2: lambda2 - lambda1 > 1
    for lambda1 in range(1, n - 1):
        for lambda2 in range(lambda1 + 2, n):
            state = np.zeros((length, 1)).astype('complex')
            k0 = 2 * np.pi * (lambda1 + lambda2) / n
            eigen = False
            guess = 1
            while (not eigen):
                func = lambda k1: (np.power(np.tan(k1 / 2), -1) -
                                   np.power(np.tan((k0 - k1) / 2), -1) -
                                   2 * np.power(np.tan(n * k1 / 2), -1))
                k1 = fsolve(func, guess)
                if k1 > 2 * np.pi:
                    k1 = k1 - 2 * np.pi
                if k1 < 0:
                    k1 = k1 + 2 * np.pi
                k2 = k0 - k1
                theta = 2 * np.pi * lambda2 - n * k2
                # print(k1,k2) # print values if needed
                count = 0
                for n1 in range(1, n):
                    for n2 in range(n1 + 1, n + 1):
                        state[count, 0] = (np.exp(1j * (k1 * n1 + k2 * n2 + 0.5 * theta)) +
                                            np.exp(1j * (k1 * n2 + k2 * n1 - 0.5 * theta)))
                        count += 1
                state = change_order_matrix(clean(normalize(state), 1e-10))
                e_state = J * (1 - np.cos(k1)) + J * (1 - np.cos(k2)) + ground_energy
                eigen = is_eigen(hamiltonian_2magnons, e_state, state)
                orth = True
                for v in states_result:
                    if (not orthogonal(v, state)):
                        orth = False
                if (not orth):
                    eigen = False
                if (eigen):
                    states_result.append(state)
                    energies.append(e_state)
                    k.append(np.array([k1, k2]))
                    # print(theta)
                elif (guess < 10):
                    guess = guess + 1
                else:
                    print("Unable to find eigenvector: ", lambda1, lambda2)
                    states_result.append(state)
                    energies.append(e_state)
                    k.append(np.array([k1, k2]))
                    break

    # class C3
    for lambda_1 in range(1, n):
        for lambda_2 in range(lambda_1, lambda_1 + 2):
            k0 = 2 * np.pi * (lambda_1 + lambda_2) / n
            phi = np.pi * (lambda_1 - lambda_2)
            eigen = False
            guess = 3
            while (not eigen):
                func = lambda nu: ((np.cos(k0 / 2)) * np.sinh(n * nu) -
                                   np.sinh((n - 1) * nu) - np.cos(phi) * np.sinh(nu))
                nu = fsolve(func, guess, xtol=1e-12, maxfev = 10000)
                # nu = anderson(func, guess)  # choose other method if needed
                # nu = newton_krylov(func, guess)  # choose other method if needed
                if abs(nu) > 0:
                    k1 = k0 / 2 + 1j * nu
                    k2 = k0 / 2 - 1j * nu
                    theta = phi + 1j * n * nu
                    vec = np.zeros((length, 1)).astype('complex')
                    count = 0
                    for n1 in range(1, n):
                        for n2 in range(n1 + 1, n + 1):
                            vec[count, 0] = (np.exp(1j * (k1 * n1 + k2 * n2 + 0.5 * theta)) +
                                             np.exp(1j * (k1 * n2 + k2 * n1 - 0.5 * theta)))
                            count = count + 1
                    vec = change_order_matrix(clean(normalize(vec), 1e-10))
                    e_state = 2*J * (1 - np.cos(k1.real) * np.cosh(k1.imag)) + ground_energy
                    eigen = is_eigen(hamiltonian_2magnons, e_state, vec)
                    orth = True
                    for v in states_result:
                        if (not orthogonal(v, vec)):
                            orth = False
                    if (not orth):
                        eigen = False
                    if (eigen):
                        states_result.append(vec)
                        energies.append(e_state)
                        k.append(np.array([k1, k2]))
                        # print(theta)
                    elif (guess < 7):
                        guess = guess + 0.2  # change initial guess
                    else:
                        # print("Unable to find eigenvector: ", lambda_1, lambda_2)
                        # print([k1, k2])
                        states_result.append(vec)
                        energies.append(e_state)
                        k.append(np.array([k1, k2]))
                        break
                else:
                    break

    if len(states_result) < int(n * (n - 1) / 2):
        vec = np.zeros((length, 1)).astype('complex')
        count = 0
        for n1 in range(1, n):
            for n2 in range(n1 + 1, n + 1):
                delta_funct = delta_function(n2, n1 + 1)
                vec[count, 0] = 1j * ((-1) ** n1) * delta_funct
                if n1 == 1 and n2 == n:
                    vec[count, 0] = 1j
                count = count + 1
        energies.append(J * (1 - np.cos((np.pi) / 2)) + ground_energy)
        states_result.append(change_order_matrix(clean(normalize(vec), 1e-10)))
        k.append(np.array([0, 0]))

    return k, states_result, energies
